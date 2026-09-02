#!/usr/bin/env python3

import math
import sys
import time


MODEL = "biohub/ESMFold2-Fast"
REVISION = "c6c7958d63f5f2f1f0fed0bb9462316f8ccceea6"
ESMC_MODEL = "biohub/ESMC-6B"
ESMC_REVISION = "650e98b7e8ed7d05a4fcf10f3b9121a55c97d23d"


def prefetch(cache_dir):
    from huggingface_hub import snapshot_download

    for repo, rev in ((MODEL, REVISION), (ESMC_MODEL, ESMC_REVISION)):
        print(f"INFO: prefetch {repo}@{rev} -> {cache_dir}", file=sys.stderr)
        snapshot_download(
            repo_id=repo,
            revision=rev,
            cache_dir=cache_dir,
            allow_patterns=["*.json", "*.safetensors"],
        )


def fasta_records(path):
    seq_id = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(sequence)
                seq_id = line[1:].split()[0]
                sequence = []
            elif seq_id is None:
                raise ValueError("sequence before first FASTA header")
            else:
                sequence.append(line)
    if seq_id is not None:
        yield seq_id, "".join(sequence)


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--prefetch":
        prefetch(sys.argv[2])
        return
    if len(sys.argv) != 4:
        raise SystemExit(
            f"usage: {sys.argv[0]} INPUT_FASTA OUTPUT_TSV THREADS\n"
            f"       {sys.argv[0]} --prefetch CACHE_DIR"
        )

    import torch
    from esm.models.esmfold2 import EsmFold2Model
    from esm.models.esmc import EsmcModel

    input_fasta, output_tsv, threads = sys.argv[1:]
    torch.set_num_threads(max(1, int(threads)))
    loaded = time.monotonic()
    model = EsmFold2Model.from_pretrained(
        MODEL,
        revision=REVISION,
        load_esmc=False,
        device="cpu",
        dtype=torch.float32,
    ).eval()
    if model.esmc is None:
        model.esmc = EsmcModel.from_pretrained(
            ESMC_MODEL,
            revision=ESMC_REVISION,
            device="cpu",
            dtype=torch.float32,
        )
    model.set_esmc_precision("fp32")
    model.set_kernel_backend(None)
    if model.device.type != "cpu":
        raise RuntimeError(f"ESMFold2-Fast loaded on {model.device}, expected CPU")
    print(f"INFO: ESMFold2-Fast model loaded in {time.monotonic() - loaded:.2f}s", file=sys.stderr)

    inference_started = time.monotonic()
    protein_count = 0
    residue_count = 0
    with open(output_tsv, "w", encoding="utf-8") as output:
        for seq_id, sequence in fasta_records(input_fasta):
            torch.manual_seed(0)
            with torch.inference_mode():
                result = model.infer_protein(
                    sequence,
                    num_loops=3,
                    num_sampling_steps=50,
                    num_diffusion_samples=1,
                )
            score = float(result["plddt"].mean()) * 100.0
            if not math.isfinite(score) or not 0.0 <= score <= 100.0:
                raise RuntimeError(f"invalid pLDDT for sequence {seq_id}: {score}")
            output.write(f"{seq_id}\t{score:.4f}\n")
            output.flush()
            protein_count += 1
            residue_count += len(sequence)
    elapsed = time.monotonic() - inference_started
    print(
        f"INFO: ESMFold2-Fast inference: {protein_count} proteins, "
        f"{residue_count} residues in {elapsed:.2f}s",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
