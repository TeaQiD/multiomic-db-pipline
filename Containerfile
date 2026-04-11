FROM python:3.12-slim

WORKDIR /app

RUN pip install --no-cache-dir \
    requests \
    pandas \
    numpy \
    scipy \
    scikit-learn \
    cptac

COPY Makefile ./
COPY scripts/ ./scripts/

CMD ["make", "-j1", "all"]
