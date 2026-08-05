""" """

from functools import cached_property
from dataclasses import dataclass
from glob import glob
from typing import Dict, Iterable, Iterator, List, Tuple

import gzip
import os
import re


from .crisprme2_error import Crisprme2AssemblyError
from .fasta_utils import FASTAEXTENSIONS, FASTA_COMPRESSED_SUFFIXES


_PANSN_SEP = "#"
_CHAIN_PREFIX_RE = re.compile(r"^(?P<sample>.+?)_hap(?P<hap>\d+)_")
_COMPRESSED = tuple(f".{s}" for s in FASTA_COMPRESSED_SUFFIXES)


@dataclass(frozen=True)
class Haplotype:
    sample_id: str
    hap_id: int  # opaque integer id from PanSN; no pat/mat meaning
    fasta: str
    chain: str


@dataclass
class AssemblySample:
    sample_id: str
    haplotypes: List[Haplotype]

    @property
    def ploidy(self) -> int:
        return len(self.haplotypes)


@dataclass
class AssemblyInputs:
    samples: List[AssemblySample]

    @cached_property
    def sample_table(self) -> "SampleTable":
        return SampleTable.from_assemblies(self)

    @property
    def ploidy(self) -> int:
        return self.samples[0].ploidy if self.samples else 0

    @property
    def sample_ids(self) -> List[str]:
        return [s.sample_id for s in self.samples]

    def haplotypes(self) -> Iterator[Haplotype]:
        for sample in self.samples:
            yield from sample.haplotypes

    @classmethod
    def discover(cls, assemblies_dir: str, chains_dir: str) -> "AssemblyInputs":
        fastas, chains = _retrieve_data(assemblies_dir, chains_dir)
        fasta_by_key = _split_fasta_by_key(fastas)
        chain_by_key = _split_chain_by_key(chains)
        _check_missing_data(fasta_by_key, chain_by_key)
        by_sample = _split_by_sample(fasta_by_key, chain_by_key)
        samples = [
            AssemblySample(sid, sorted(haps, key=lambda h: h.hap_id))
            for sid, haps in sorted(by_sample.items())
        ]
        ploidies = {s.ploidy for s in samples}
        if len(ploidies) > 1:
            detail = ", ".join(f"{s.sample_id}={s.ploidy}" for s in samples)
            raise Crisprme2AssemblyError(
                f"Non-uniform ploidy across samples ({detail}); every sample in "
                "a run must have the same number of haplotypes"
            )
        return cls(samples)


class SampleTable:

    def __init__(self, names: Iterable[str]) -> None:
        self._names: Tuple[str, ...] = tuple(sorted(dict.fromkeys(names)))
        self._index: Dict[str, int] = {n: i for i, n in enumerate(self._names)}

    @classmethod
    def from_assemblies(cls, assemblies: "AssemblyInputs") -> "SampleTable":
        return cls(assemblies.sample_ids)

    def index(self, sample_id: str) -> int:
        try:
            return self._index[sample_id]
        except KeyError:
            raise Crisprme2AssemblyError(f"Unknown sample id {sample_id!r}")

    def name(self, index: int) -> str:
        try:
            return self._names[index]
        except IndexError:
            raise Crisprme2AssemblyError(f"Sample index {index} out of range")

    @property
    def names(self) -> Tuple[str, ...]:
        """Index → name, in u32-index order (what the native table consumes)."""
        return self._names

    def __len__(self) -> int:
        return len(self._names)

    def __iter__(self) -> Iterator[str]:
        return iter(self._names)


def _glob_fastas(folder: str) -> List[str]:
    files: List[str] = []
    for ext in FASTAEXTENSIONS:
        files.extend(glob(os.path.join(folder, f"*.{ext}")))
        for suf in FASTA_COMPRESSED_SUFFIXES:
            files.extend(glob(os.path.join(folder, f"*.{ext}.{suf}")))
    return sorted(set(files))


def _glob_chains(folder: str) -> List[str]:
    files: List[str] = []
    for pattern in ("*.chain", "*.chain.gz"):
        files.extend(glob(os.path.join(folder, pattern)))
    return sorted(set(files))


def _retrieve_data(assemblies_dir: str, chains_dir: str) -> Tuple[List[str], List[str]]:
    fastas = _glob_fastas(assemblies_dir)
    if not fastas:
        raise Crisprme2AssemblyError(f"No FASTA file found in {assemblies_dir}")
    chains = _glob_chains(chains_dir)
    if not chains:
        raise Crisprme2AssemblyError(f"No chain file found in {chains_dir}")
    return fastas, chains


def _parse_pansn_header(header: str) -> Tuple[str, int]:
    parts = header.split(_PANSN_SEP, 2)  # tolerate '#' inside the contig field
    if len(parts) != 3 or not parts[0] or not parts[1].isdigit():
        raise Crisprme2AssemblyError(
            f"Non-conforming FASTA header {header!r}: expected PanSN "
            "'sample#haplotype#contig' with an integer haplotype id"
        )
    return parts[0], int(parts[1])


def _first_fasta_header(path: str) -> str:
    opener = gzip.open if path.endswith(_COMPRESSED) else open
    try:
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    return line[1:].strip().split()[0]
    except OSError as e:
        raise Crisprme2AssemblyError(f"Cannot read FASTA {path}: {e}") from e
    raise Crisprme2AssemblyError(f"No FASTA header found in {path}")


def _split_fasta_by_key(fastas: List[str]) -> Dict[Tuple[str, int], str]:
    fasta_by_key: Dict[Tuple[str, int], str] = {}
    for fa in fastas:
        pansn_key = _parse_pansn_header(_first_fasta_header(fa))
        if pansn_key in fasta_by_key:
            raise Crisprme2AssemblyError(
                f"Duplicate haplotype {pansn_key[0]}#hap{pansn_key[1]} in "
                f"assemblies: {fasta_by_key[pansn_key]} and {fa}"
            )
        fasta_by_key[pansn_key] = fa
    return fasta_by_key


def _parse_chain_filename(path: str) -> Tuple[str, int]:
    name = os.path.basename(path)
    m = _CHAIN_PREFIX_RE.match(name)
    if not m:
        raise Crisprme2AssemblyError(
            f"Chain filename {name!r} does not match the required prefix "
            "'<sample-id>_hap<hap-id>_'"
        )
    return m.group("sample"), int(m.group("hap"))


def _split_chain_by_key(chains: List[str]) -> Dict[Tuple[str, int], str]:
    chain_by_key: Dict[Tuple[str, int], str] = {}
    for ch in chains:
        pansn_key = _parse_chain_filename(ch)
        if pansn_key in chain_by_key:
            raise Crisprme2AssemblyError(
                f"Duplicate haplotype {pansn_key[0]}#hap{pansn_key[1]} in "
                f"chains: {chain_by_key[pansn_key]} and {ch}"
            )
        chain_by_key[pansn_key] = ch
    return chain_by_key


def _check_missing_data(
    fasta_by_key: Dict[Tuple[str, int], str], chain_by_key: Dict[Tuple[str, int], str]
) -> None:
    missing_chain = sorted(fasta_by_key.keys() - chain_by_key.keys())
    missing_fasta = sorted(chain_by_key.keys() - fasta_by_key.keys())
    if missing_chain:
        raise Crisprme2AssemblyError(
            "Missing chain file for haplotype(s): "
            + ", ".join(f"{s}#hap{h}" for s, h in missing_chain)
        )
    if missing_fasta:
        raise Crisprme2AssemblyError(
            "Chain file(s) without a matching assembly FASTA: "
            + ", ".join(f"{s}#hap{h}" for s, h in missing_fasta)
        )


def _split_by_sample(
    fasta_by_key: Dict[Tuple[str, int], str], chain_by_key: Dict[Tuple[str, int], str]
) -> Dict[str, List[Haplotype]]:
    by_sample: Dict[str, List[Haplotype]] = {}
    for (sample_id, hap_id), fa in fasta_by_key.items():
        by_sample.setdefault(sample_id, []).append(
            Haplotype(sample_id, hap_id, fa, chain_by_key[(sample_id, hap_id)])
        )
    return by_sample


def _haplotype_key_of_contigs(contigs: List[str], fasta: str) -> Tuple[str, int]:
    if not contigs:
        raise Crisprme2AssemblyError(f"FASTA {fasta} has no contigs")
    keys = {_parse_pansn_header(c) for c in contigs}
    if len(keys) != 1:
        offenders = ", ".join(sorted(f"{s}#hap{h}" for s, h in keys))
        raise Crisprme2AssemblyError(
            f"FASTA {fasta} mixes multiple haplotypes ({offenders}); one "
            "haplotype (sample#hap) per file is required"
        )
    return next(iter(keys))


def validate_haplotype_contigs(haplotype: Haplotype, contigs: List[str]) -> None:
    key = _haplotype_key_of_contigs(contigs, haplotype.fasta)
    if key != (haplotype.sample_id, haplotype.hap_id):
        raise Crisprme2AssemblyError(
            f"FASTA {haplotype.fasta} headers report {key[0]}#hap{key[1]} but "
            f"was registered as {haplotype.sample_id}#hap{haplotype.hap_id}"
        )
