default: test

# Fast test mode (~5-10 min, ~200 MB) — limits TCGA fetches, keeps full CPTAC
test:
    PIPELINE_MODE=test podman compose up --build --abort-on-container-exit
    podman compose down

# Full pipeline (~3-4 hr, ~12 GB DB) — fetches every TCGA file
full:
    PIPELINE_MODE=full podman compose up --build --abort-on-container-exit
    podman compose down

clean:
    podman compose down --rmi all 2>/dev/null || true
    rm -f data/*.db data/*.csv data/.*.done
    rm -f figures/*.png
