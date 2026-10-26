"""Release smoke tests; no Docker daemon, builds, or registry access required."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class DockerRelease(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="docker release ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.context = self.root / "docker"
        shutil.copytree(Path(__file__).resolve().parent, self.context,
                        ignore=shutil.ignore_patterns("__pycache__"))
        self.log = self.root / "calls.jsonl"
        # Keep executable and build contexts in separate directories.
        executable = self.root / "bin" / "docker"
        executable.parent.mkdir()
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import json, os, pathlib, sys\n"
            "log = pathlib.Path(os.environ['DOCKER_TEST_LOG'])\n"
            "with log.open('a') as handle:\n"
            "    handle.write(json.dumps(sys.argv[1:]) + '\\n')\n"
            "count = len(log.read_text().splitlines())\n"
            "sys.exit(42 if str(count) == os.environ.get('FAIL_CALL') else 0)\n"
        )
        executable.chmod(0o755)
        self.env = dict(os.environ, PATH=str(executable.parent) + os.pathsep + os.environ['PATH'],
                        DOCKER_TEST_LOG=str(self.log))

    def run_release(self, args, fail_call=None):
        env = dict(self.env)
        if fail_call is not None:
            env['FAIL_CALL'] = str(fail_call)
        result = subprocess.run(['sh', str(self.context / 'update_docker.sh'), *args],
                                cwd=self.root, env=env, capture_output=True, text=True)
        calls = [json.loads(line) for line in self.log.read_text().splitlines()] if self.log.exists() else []
        return result, calls

    def test_release_version_and_tags_from_unrelated_directory(self):
        result, calls = self.run_release(['0.9.0-rc1'])
        self.assertEqual(result.returncode, 0, result.stderr)
        env_image = 'naechyun/leviosam2_env:0.9.0-rc1'
        self.assertEqual(calls, [
            ['build', '--build-arg', 'VERSION=0.9.0-rc1', '--tag', env_image, str(self.context / 'env')],
            ['build', '--build-arg', 'ENV_IMAGE=' + env_image,
             '--tag', 'naechyun/leviosam2:0.9.0-rc1', '--tag', 'naechyun/leviosam2:latest',
             str(self.context / 'production')],
            ['push', 'naechyun/leviosam2:0.9.0-rc1'],
            ['push', 'naechyun/leviosam2:latest'],
        ])
        dockerfile = (self.context / 'env' / 'Dockerfile').read_text()
        self.assertIn('ARG VERSION\n', dockerfile)
        self.assertIn('refs/tags/v${VERSION}.tar.gz', dockerfile)
        self.assertNotIn('ENV VERSION', dockerfile)
        self.assertIn('FROM ${ENV_IMAGE} AS lev',
                      (self.context / 'production' / 'Dockerfile').read_text())

    def test_build_and_push_failures_stop_immediately(self):
        for failure in (1, 2, 3, 4):
            with self.subTest(failure=failure):
                self.log.unlink(missing_ok=True)
                result, calls = self.run_release(['0.9.0'], failure)
                self.assertEqual(result.returncode, 42)
                self.assertEqual(len(calls), failure)
                if failure <= 2:
                    self.assertTrue(all(call[0] == 'build' for call in calls))

    def test_invalid_arguments_do_not_build_or_push(self):
        for args in ([], [''], ['v0.9.0'], ['0.9/invalid'], ['0.9.0', 'extra']):
            with self.subTest(args=args):
                result, calls = self.run_release(args)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(calls, [])


if __name__ == '__main__':
    unittest.main()
