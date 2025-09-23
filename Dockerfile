# ---- Base image ----
FROM python:3.11-slim-bullseye AS base

# Prevent Python from writing .pyc files; unbuffer logs
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Faster, repeatable installs
ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_PROGRESS_BAR=off

# System deps: tini (init), certs, ps (procps), and bash for scripts
RUN apt-get update && apt-get install -y --no-install-recommends \
      tini ca-certificates procps bash \
    && rm -rf /var/lib/apt/lists/*

# Non-root user
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} app && useradd -m -u ${UID} -g ${GID} app
WORKDIR /app

# ---- Builder layer (for caching installs) ----
FROM base AS builder

# Only copy files that affect dependency resolution first (better caching)
COPY pyproject.toml README.md LICENSE /app/
# If you have a src layout, copy it next
COPY src /app/src

# Install the package (and its runtime deps) into a wheels dir
RUN python -m pip install --upgrade pip wheel \
 && python -m pip wheel --no-deps --wheel-dir /wheels /app

# ---- Final runtime image ----
FROM base AS runtime

# Copy built wheels and install
COPY --from=builder /wheels /wheels
RUN python -m pip install /wheels/*.whl

# Drop to non-root
USER app

# Default workdir where users can mount data
WORKDIR /workspace

# Use tini as minimal init to handle signals properly
ENTRYPOINT ["/usr/bin/tini", "--"]

# Expose the CLI by default; override with `docker run ... polycore --help`
CMD ["polycore", "--help"]
