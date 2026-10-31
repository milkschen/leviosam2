"""Restore original library identity after FASTQ realignment, using samtools.

Only Python's standard library is required. The read-name index is kept on disk
so memory use does not grow with the number of deferred reads.
"""

import argparse
from contextlib import closing
import itertools
import pathlib
import shlex
import sqlite3
import subprocess
import sys
import tempfile


def read_header(stream):
    header = []
    for line in stream:
        if line.startswith('@'):
            header.append(line.rstrip('\n'))
        else:
            return header, itertools.chain([line], stream)
    return header, iter(())


def read_groups(header):
    groups = {}
    for line in header:
        if not line.startswith('@RG\t'):
            continue
        fields = {}
        for field in line.split('\t')[1:]:
            key, separator, value = field.partition(':')
            if not separator or key in fields:
                raise ValueError(f'Invalid read-group header: {line}')
            fields[key] = value
        identity = fields.get('ID')
        if not identity or identity in groups:
            raise ValueError(f'Missing or duplicate read-group ID: {line}')
        groups[identity] = fields
    return groups


def merge_groups(*collections):
    groups = {}
    for collection in collections:
        for identity, fields in collection.items():
            if identity in groups and groups[identity] != fields:
                raise ValueError(f'Conflicting definitions for read-group ID {identity!r}')
            groups[identity] = fields
    return groups


def sam_record(line, paired):
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 11:
        raise ValueError('Malformed SAM record while restoring read groups')
    flag = int(fields[1])
    end = flag & (0x40 | 0x80) if paired else 0
    if paired and end not in (0x40, 0x80):
        raise ValueError(f'Paired record {fields[0]!r} has no unique read end')
    tags = [field for field in fields[11:] if field.startswith('RG:')]
    if len(tags) > 1 or (tags and not tags[0].startswith('RG:Z:')):
        raise ValueError(f'Invalid or duplicate RG tag on {fields[0]!r}')
    group = tags[0][5:] if tags else None
    return fields, flag, end, group


def samtools_header(samtools, path):
    return subprocess.check_output(
        [*samtools, 'view', '-H', str(path)], text=True).splitlines()


def index_source(database, samtools, source, paired):
    database.execute('CREATE TABLE reads (name TEXT, end INTEGER, rg TEXT, '
                     'PRIMARY KEY (name, end)) WITHOUT ROWID')
    process = subprocess.Popen([*samtools, 'view', '-h', str(source)],
                               stdout=subprocess.PIPE, text=True)
    try:
        header, lines = read_header(process.stdout)
        groups = read_groups(header)
        for line in lines:
            fields, flag, end, group = sam_record(line, paired)
            if flag & (0x100 | 0x800):
                continue
            if group is not None and group not in groups:
                raise ValueError(f'Undefined original read-group {group!r} on {fields[0]!r}')
            key = (fields[0], end)
            cursor = database.execute('INSERT OR IGNORE INTO reads VALUES (?, ?, ?)',
                                      (*key, group))
            if not cursor.rowcount:
                existing = database.execute('SELECT rg FROM reads WHERE name=? AND end=?', key).fetchone()
                if existing[0] != group:
                    raise ValueError(f'Ambiguous original read-group for read {fields[0]!r}, end {end}')
        if process.wait():
            raise RuntimeError(f'samtools could not read original read groups from {source}')
        database.commit()
        return groups
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.kill()
            process.wait()


def restore(samtools, source, paired, stream, output, header_from=()):
    with tempfile.TemporaryDirectory(prefix='leviosam2-read-groups-') as directory:
        with closing(sqlite3.connect(str(pathlib.Path(directory) / 'reads.sqlite'))) as database:
            database.execute('PRAGMA journal_mode=OFF')
            database.execute('PRAGMA synchronous=OFF')
            source_groups = index_source(database, samtools, source, paired)
            header, lines = read_header(stream)
            groups = merge_groups(read_groups(header), source_groups,
                *(read_groups(samtools_header(samtools, path)) for path in header_from))
            # Preserve target sequence/program headers; add all library definitions.
            for line in header:
                if not line.startswith('@RG\t'):
                    output.write(line + '\n')
            for fields in groups.values():
                output.write('@RG\t' + '\t'.join(f'{key}:{value}' for key, value in fields.items()) + '\n')
            for line in lines:
                fields, _, end, fallback = sam_record(line, paired)
                original = database.execute('SELECT rg FROM reads WHERE name=? AND end=?',
                                            (fields[0], end)).fetchone()
                if original is None:
                    raise ValueError(f'Realigned read {fields[0]!r}, end {end} has no source record')
                # An explicit aligner RG is a fallback for originally ungrouped reads.
                group = original[0] if original[0] is not None else fallback
                if group is not None and group not in groups:
                    raise ValueError(f'Undefined realigned read-group {group!r} on {fields[0]!r}')
                fields = fields[:11] + [field for field in fields[11:] if not field.startswith('RG:')]
                if group is not None:
                    fields.append('RG:Z:' + group)
                output.write('\t'.join(fields) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samtools', required=True)
    parser.add_argument('--source', type=pathlib.Path)
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--header-from', type=pathlib.Path, action='append', default=[])
    parser.add_argument('--check-headers', type=pathlib.Path, nargs='+')
    args = parser.parse_args()
    samtools = [args.samtools] if pathlib.Path(args.samtools).is_file() else shlex.split(args.samtools)
    if args.check_headers:
        merge_groups(*(read_groups(samtools_header(samtools, path)) for path in args.check_headers))
    elif args.source is None:
        parser.error('--source is required when restoring read groups')
    else:
        restore(samtools, args.source, args.paired, sys.stdin, sys.stdout, args.header_from)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, sqlite3.Error, subprocess.SubprocessError) as error:
        print(f'[E::restore_read_groups] {error}', file=sys.stderr)
        sys.exit(1)
