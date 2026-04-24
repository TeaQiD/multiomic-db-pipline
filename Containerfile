FROM python:3.12-slim

WORKDIR /app

# Build tools are needed because a few cptac transitive deps (ncls,
# sorted_nearest via pyranges) no longer ship wheels for python:3.12-slim
# and must compile from source.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    requests \
    pandas \
    numpy \
    scipy \
    scikit-learn \
    cptac \
    matplotlib==3.9.2 \
    seaborn==0.13.2

COPY Makefile ./
COPY scripts/ ./scripts/

CMD ["make", "-j1", "all"]
