Run `sh /path/to/leviosam2/docker/update_docker.sh 0.9.0` from any directory,
using the release version without the leading `v`. The corresponding GitHub tag
(`v0.9.0` in this example) must exist. Docker registry credentials must already be
configured for `naechyun/leviosam2`.

The script builds a versioned environment image with `--build-arg VERSION=0.9.0`,
then passes that exact image to the production build as `ENV_IMAGE`. It builds
and pushes `naechyun/leviosam2:0.9.0` and `naechyun/leviosam2:latest`, stopping
immediately if any build or push fails. The environment image stays local.

Run `python3 docker/test_update_docker.py` for smoke tests with a fake Docker
executable; these tests do not build images or push to a registry.
