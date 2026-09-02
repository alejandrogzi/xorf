# ==============================================================================
# blast.Dockerfile - Container for 'orf blast' subcommand
# ==============================================================================
# Includes: orf binary + diamond + orfipy + CPU ESMFold2-Fast

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
    git \
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

# Create Python 3.12 environment with protein-search tools
RUN micromamba create -y -n blastenv -c conda-forge -c bioconda \
    python=3.12 \
    pip \
    diamond=2.2.4 \
    mmseqs2 \
    && micromamba clean -a -y

# Install orfipy via pip (requires gcc for Cython compilation)
RUN micromamba run -n blastenv pip install --no-cache-dir \
    "orfipy>=0.0.4"

# Install PSAURON and ESM in the same transaction so one CPU-only PyTorch owns
# the environment. Upstream ESM currently declares CUDA cuequivariance wheels
# unconditionally on x86_64; remove them because ESMFold2 has a PyTorch fallback.
RUN micromamba run -n blastenv pip install --no-cache-dir \
    --extra-index-url https://download.pytorch.org/whl/cpu \
    "numpy<2" \
    "torch==2.11.0+cpu" \
    "torchvision==0.26.0+cpu" \
    "torchaudio==2.11.0+cpu" \
    "psauron==1.1.3" \
    "transformers==4.57.6" \
    "esm @ git+https://github.com/evolutionaryscale/esm.git@827ec128e4cdaf80f7d6f95fb367a08980b34918" && \
    micromamba run -n blastenv pip uninstall -y \
    cuequivariance-torch \
    cuequivariance-ops-torch-cu13 \
    cuequivariance-ops-cu13 \
    nvidia-cublas \
    nvidia-cuda-nvrtc \
    nvidia-ml-py

# Copy Rust binary
USER root
COPY --from=builder /build/target/release/orf /usr/local/bin/orf
RUN chmod +x /usr/local/bin/orf

# Set environment variables so conda environment is available
ENV PATH="/opt/conda/envs/blastenv/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/blastenv/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH="/opt/conda/envs/blastenv/lib/python3.12/site-packages:$PYTHONPATH"

# Create non-root user
RUN useradd -m -u 1000 orfuser
USER orfuser
WORKDIR /data

# Verify installations (tools should be directly in PATH now)
RUN orf blast --help && \
    diamond --version && \
    mmseqs version && \
    orfipy --version && \
    python -c "import torch, torchaudio, torchvision; from esm.models.esmfold2 import EsmFold2Model; assert torch.version.cuda is None"

# ENTRYPOINT ["orf"]
# CMD ["blast", "--help"]
