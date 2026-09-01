# ==============================================================================
# blast.Dockerfile - Container for 'orf blast' subcommand
# ==============================================================================
# Includes: orf binary + diamond + orfipy (Python 3.9)
# Result: ~500MB container

# ---------- Build Stage ----------
FROM rust:1.93-slim-bullseye AS builder

WORKDIR /build

COPY modules/orf/Cargo.toml modules/orf/Cargo.lock ./
COPY modules/orf/src/ ./src/

RUN cargo build --release && \
    strip target/release/orf

# ---------- Runtime Stage ----------
FROM mambaorg/micromamba:1.5.8

USER root

# Install system dependencies including gcc for building Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Create necessary directories with proper permissions
RUN mkdir -p /opt/conda && \
    chmod -R 0777 /opt/conda

USER $MAMBA_USER

# Set environment variables
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Create Python 3.9 environment with diamond and orfipy
RUN micromamba create -y -n blastenv -c conda-forge -c bioconda \
    python=3.9 \
    pip \
    diamond=2.2.4 \
    psauron=1.1.3 \
    mmseqs2 \
    && micromamba clean -a -y

# Install orfipy via pip (requires gcc for Cython compilation)
RUN micromamba run -n blastenv pip install --no-cache-dir \
    "orfipy>=0.0.4"

# Copy Rust binary
USER root
COPY --from=builder /build/target/release/orf /usr/local/bin/orf
RUN chmod +x /usr/local/bin/orf

# Set environment variables so conda environment is available
ENV PATH="/opt/conda/envs/blastenv/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/blastenv/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH="/opt/conda/envs/blastenv/lib/python3.9/site-packages:$PYTHONPATH"

# Create non-root user
RUN useradd -m -u 1000 orfuser
USER orfuser
WORKDIR /data

# Verify installations (tools should be directly in PATH now)
RUN orf blast --help && \
    diamond --version && \
    mmseqs version && \
    orfipy --version

# ENTRYPOINT ["orf"]
# CMD ["blast", "--help"]
