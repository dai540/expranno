# expranno Docker helper

This directory contains a minimal Docker recipe for rebuilding an
`expranno` environment with optional backends.

Build:

```bash
docker build -t expranno-local -f inst/docker/Dockerfile .
```

Run:

```bash
docker run --rm -it -v "$PWD":/work expranno-local
```

The image installs `expranno`, then calls the bundled
`expranno_install_optional_backends()` helper so the container can also
exercise hybrid annotation, deconvolution, and signature scoring
examples.
