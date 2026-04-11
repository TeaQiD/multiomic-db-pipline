default: build run

build:
    podman compose build

run:
    podman compose up --abort-on-container-exit
    podman compose down

clean:
    podman compose down --rmi all 2>/dev/null || true
    rm -f data/*.db
