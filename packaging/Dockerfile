# ZFN WebSuite — reproducible build
#
#   docker build -t zfn-websuite:0.7.0 .
#   docker run --rm -p 8797:8797 -v "$PWD/out:/out" zfn-websuite:0.7.0
#
# The image is deliberately small: the package itself is standard-library only.
# numpy and pyahocorasick are installed because they speed the whole-genome
# filters up by an order of magnitude, but the code falls back to pure Python
# when they are missing.

FROM python:3.12-slim

LABEL org.opencontainers.image.title="ZFN WebSuite" \
      org.opencontainers.image.description="Repeat-aware, genome-scale ZFN design with order-ready CDS" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="0.7.0"

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONIOENCODING=utf-8

WORKDIR /app

# 1) optional accelerators first (they are pure wheels, no compiler needed)
RUN pip install --no-cache-dir "numpy>=1.26" "pyahocorasick>=2.0"

# 2) the package itself
COPY pyproject.toml README.md ./
COPY zfn_design ./zfn_design
RUN pip install --no-cache-dir --no-build-isolation .

# 3) a non-root user, because the server binds a port and writes result dirs
RUN useradd --create-home --uid 1000 runner && mkdir -p /out && chown -R runner /out
USER runner

EXPOSE 8797
VOLUME ["/out"]

# self-check at build time: the CLI must answer without any arguments
RUN zfn-design info >/dev/null

CMD ["zfn-design", "gui", "--host", "0.0.0.0", "--port", "8797", "--no-browser", "--out-root", "/out"]
