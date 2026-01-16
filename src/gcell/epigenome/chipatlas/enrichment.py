"""ChIP-Atlas Enrichment Analysis API interface.

This module provides programmatic access to ChIP-Atlas enrichment analysis,
allowing you to find transcription factors or histone modifications
enriched in your regions of interest.
"""

from __future__ import annotations

import io
import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Sequence

# ChIP-Atlas enrichment analysis endpoint
ENRICHMENT_ENDPOINT = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"

# Supported antigen classes for enrichment analysis
ENRICHMENT_ANTIGEN_CLASSES = [
    "TFs and others",
    "Histone",
    "Input control",
    "RNA polymerase",
    "ATAC-Seq",
    "DNase-seq",
    "Bisulfite-Seq",
]

# Supported thresholds
ENRICHMENT_THRESHOLDS = [50, 100, 200, 500]


@dataclass
class EnrichmentResult:
    """Result from ChIP-Atlas enrichment analysis.

    Attributes
    ----------
    request_id : str
        Unique identifier for this analysis request
    status : str
        Current status of the analysis
    results : pd.DataFrame | None
        Results DataFrame when analysis is complete
    error : str | None
        Error message if analysis failed
    """

    request_id: str
    status: str = "pending"
    results: pd.DataFrame | None = None
    error: str | None = None
    _raw_response: dict = field(default_factory=dict, repr=False)

    @property
    def is_complete(self) -> bool:
        """Check if analysis is complete."""
        return self.status in ("finished", "error")

    @property
    def is_successful(self) -> bool:
        """Check if analysis completed successfully."""
        return self.status == "finished" and self.results is not None


class ChipAtlasEnrichment:
    """Interface for ChIP-Atlas enrichment analysis.

    This class provides programmatic access to ChIP-Atlas enrichment analysis,
    which identifies transcription factors, histone modifications, or other
    epigenetic features enriched in your genomic regions.

    Parameters
    ----------
    timeout : int, optional
        Request timeout in seconds. Default is 30.
    poll_interval : int, optional
        Interval between status checks in seconds. Default is 5.

    Examples
    --------
    >>> from gcell.epigenome.chipatlas import ChipAtlasEnrichment
    >>> enrichment = ChipAtlasEnrichment()

    Run enrichment analysis with BED data:

    >>> # Define your regions of interest
    >>> bed_data = '''chr1\\t1000\\t2000
    ... chr1\\t3000\\t4000
    ... chr2\\t5000\\t6000'''
    >>> result = enrichment.analyze_regions(
    ...     bed_data,
    ...     assembly="hg38",
    ...     antigen_class="TFs and others",
    ... )

    Run enrichment with gene list:

    >>> genes = ["TP53", "MYC", "BRCA1", "PTEN"]
    >>> result = enrichment.analyze_genes(
    ...     genes,
    ...     assembly="hg38",
    ...     distance_up=5000,
    ...     distance_down=5000,
    ... )
    """

    def __init__(
        self,
        timeout: int = 30,
        poll_interval: int = 5,
    ):
        self.timeout = timeout
        self.poll_interval = poll_interval

    def analyze_regions(
        self,
        bed_data: str | Path | pd.DataFrame,
        assembly: str,
        antigen_class: str = "TFs and others",
        cell_class: str = "All cell types",
        threshold: int = 200,
        compare_with: Literal["random", "bed"] = "random",
        compare_bed: str | Path | pd.DataFrame | None = None,
        permutations: int = 1,
        title: str = "",
        wait: bool = True,
        max_wait: int = 600,
    ) -> EnrichmentResult:
        """Run enrichment analysis on genomic regions.

        Parameters
        ----------
        bed_data : str | Path | pd.DataFrame
            BED format data as string, file path, or DataFrame
        assembly : str
            Genome assembly (hg38, hg19, mm10, etc.)
        antigen_class : str
            Antigen class to analyze
        cell_class : str
            Cell type class to analyze
        threshold : int
            Peak significance threshold (50, 100, 200, or 500)
        compare_with : str
            Comparison method: "random" or "bed"
        compare_bed : str | Path | pd.DataFrame, optional
            BED data to compare against (required if compare_with="bed")
        permutations : int
            Number of permutations (1, 10, or 100)
        title : str, optional
            Analysis title
        wait : bool
            If True, wait for analysis to complete
        max_wait : int
            Maximum time to wait in seconds

        Returns
        -------
        EnrichmentResult
            Analysis result object
        """
        # Prepare BED data
        bed_str = self._prepare_bed_data(bed_data)

        # Prepare comparison data
        if compare_with == "random":
            type_b = "rnd"
            bed_b = "empty"
        elif compare_with == "bed":
            if compare_bed is None:
                raise ValueError("compare_bed is required when compare_with='bed'")
            type_b = "bed"
            bed_b = self._prepare_bed_data(compare_bed)
        else:
            raise ValueError(f"Invalid compare_with: {compare_with}")

        # Validate parameters
        if threshold not in ENRICHMENT_THRESHOLDS:
            raise ValueError(
                f"Invalid threshold: {threshold}. "
                f"Must be one of {ENRICHMENT_THRESHOLDS}"
            )

        if permutations not in [1, 10, 100]:
            raise ValueError("permutations must be 1, 10, or 100")

        # Build request
        params = {
            "format": "text",
            "result": "www",
            "genome": assembly,
            "antigenClass": antigen_class,
            "cellClass": cell_class,
            "threshold": str(threshold),
            "typeA": "bed",
            "bedAFile": bed_str,
            "typeB": type_b,
            "bedBFile": bed_b,
            "permTime": str(permutations),
        }

        if title:
            params["title"] = title

        # Submit request
        result = self._submit_request(params)

        if wait and not result.is_complete:
            result = self._wait_for_completion(result, max_wait)

        return result

    def analyze_genes(
        self,
        genes: Sequence[str],
        assembly: str,
        antigen_class: str = "TFs and others",
        cell_class: str = "All cell types",
        threshold: int = 200,
        distance_up: int = 5000,
        distance_down: int = 5000,
        compare_with: Literal["random", "refseq", "genes"] = "random",
        compare_genes: Sequence[str] | None = None,
        permutations: int = 1,
        title: str = "",
        wait: bool = True,
        max_wait: int = 600,
    ) -> EnrichmentResult:
        """Run enrichment analysis on gene list.

        Parameters
        ----------
        genes : Sequence[str]
            List of gene symbols
        assembly : str
            Genome assembly (hg38, hg19, mm10, etc.)
        antigen_class : str
            Antigen class to analyze
        cell_class : str
            Cell type class to analyze
        threshold : int
            Peak significance threshold (50, 100, 200, or 500)
        distance_up : int
            Distance upstream of TSS in bp
        distance_down : int
            Distance downstream of TSS in bp
        compare_with : str
            Comparison method: "random", "refseq", or "genes"
        compare_genes : Sequence[str], optional
            Gene list to compare against (required if compare_with="genes")
        permutations : int
            Number of permutations (1, 10, or 100)
        title : str, optional
            Analysis title
        wait : bool
            If True, wait for analysis to complete
        max_wait : int
            Maximum time to wait in seconds

        Returns
        -------
        EnrichmentResult
            Analysis result object
        """
        # Prepare gene list
        gene_str = "\n".join(genes)

        # Prepare comparison
        if compare_with == "random":
            type_b = "rnd"
            bed_b = "empty"
        elif compare_with == "refseq":
            type_b = "refseq"
            bed_b = "empty"
        elif compare_with == "genes":
            if compare_genes is None:
                raise ValueError("compare_genes is required when compare_with='genes'")
            type_b = "userlist"
            bed_b = "\n".join(compare_genes)
        else:
            raise ValueError(f"Invalid compare_with: {compare_with}")

        # Validate parameters
        if threshold not in ENRICHMENT_THRESHOLDS:
            raise ValueError(
                f"Invalid threshold: {threshold}. "
                f"Must be one of {ENRICHMENT_THRESHOLDS}"
            )

        # Build request
        params = {
            "format": "text",
            "result": "www",
            "genome": assembly,
            "antigenClass": antigen_class,
            "cellClass": cell_class,
            "threshold": str(threshold),
            "typeA": "gene",
            "bedAFile": gene_str,
            "typeB": type_b,
            "bedBFile": bed_b,
            "permTime": str(permutations),
            "distanceUp": str(distance_up),
            "distanceDown": str(distance_down),
        }

        if title:
            params["title"] = title

        # Submit request
        result = self._submit_request(params)

        if wait and not result.is_complete:
            result = self._wait_for_completion(result, max_wait)

        return result

    def _prepare_bed_data(self, data: str | Path | pd.DataFrame) -> str:
        """Convert BED data to string format."""
        if isinstance(data, pd.DataFrame):
            # Assume standard BED columns
            cols = data.columns[:3]  # chrom, start, end
            lines = []
            for _, row in data.iterrows():
                lines.append(f"{row[cols[0]]}\t{row[cols[1]]}\t{row[cols[2]]}")
            return "\n".join(lines)
        elif isinstance(data, Path) or (isinstance(data, str) and Path(data).exists()):
            return Path(data).read_text()
        else:
            # Assume it's already a string
            return str(data).replace(" ", "_")

    def _submit_request(self, params: dict) -> EnrichmentResult:
        """Submit enrichment analysis request."""
        data = urllib.parse.urlencode(params).encode("utf-8")

        req = urllib.request.Request(
            ENRICHMENT_ENDPOINT,
            data=data,
            method="POST",
        )

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                response_text = response.read().decode("utf-8")
        except Exception as e:
            return EnrichmentResult(
                request_id="",
                status="error",
                error=f"Failed to submit request: {e}",
            )

        # Parse response - ChIP-Atlas returns the request ID
        request_id = response_text.strip()

        return EnrichmentResult(
            request_id=request_id,
            status="submitted",
        )

    def check_status(self, result: EnrichmentResult) -> EnrichmentResult:
        """Check the status of an enrichment analysis.

        Parameters
        ----------
        result : EnrichmentResult
            Result object from a previous analysis

        Returns
        -------
        EnrichmentResult
            Updated result object
        """
        if not result.request_id:
            return result

        url = f"{ENRICHMENT_ENDPOINT}{result.request_id}"

        try:
            with urllib.request.urlopen(url, timeout=self.timeout) as response:
                response_text = response.read().decode("utf-8")
        except Exception as e:
            result.error = f"Failed to check status: {e}"
            return result

        # Parse status response
        try:
            status_data = json.loads(response_text)
            result._raw_response = status_data
            result.status = status_data.get("status", "unknown")

            if result.status == "finished":
                # Fetch results
                result = self._fetch_results(result)

        except json.JSONDecodeError:
            # May be plain text status
            if "finished" in response_text.lower():
                result.status = "finished"
                result = self._fetch_results(result)
            elif "running" in response_text.lower():
                result.status = "running"
            elif "error" in response_text.lower():
                result.status = "error"
                result.error = response_text

        return result

    def _fetch_results(self, result: EnrichmentResult) -> EnrichmentResult:
        """Fetch results for a completed analysis."""
        url = f"{ENRICHMENT_ENDPOINT}{result.request_id}?info=result&format=tsv"

        try:
            with urllib.request.urlopen(url, timeout=self.timeout) as response:
                content = response.read().decode("utf-8")

            result.results = pd.read_csv(io.StringIO(content), sep="\t")
            result.status = "finished"

        except Exception as e:
            result.error = f"Failed to fetch results: {e}"

        return result

    def _wait_for_completion(
        self,
        result: EnrichmentResult,
        max_wait: int,
    ) -> EnrichmentResult:
        """Wait for analysis to complete."""
        start_time = time.time()

        while not result.is_complete:
            if time.time() - start_time > max_wait:
                result.status = "timeout"
                result.error = f"Analysis timed out after {max_wait} seconds"
                break

            time.sleep(self.poll_interval)
            result = self.check_status(result)

        return result

    def get_result(
        self,
        request_id: str,
        wait: bool = True,
        max_wait: int = 600,
    ) -> EnrichmentResult:
        """Get results for a previously submitted analysis.

        Parameters
        ----------
        request_id : str
            Request ID from a previous submission
        wait : bool
            If True, wait for analysis to complete
        max_wait : int
            Maximum time to wait in seconds

        Returns
        -------
        EnrichmentResult
            Analysis result object
        """
        result = EnrichmentResult(request_id=request_id, status="checking")
        result = self.check_status(result)

        if wait and not result.is_complete:
            result = self._wait_for_completion(result, max_wait)

        return result
