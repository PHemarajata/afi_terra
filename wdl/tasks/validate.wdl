version 1.0

task CompareExpectedConcordance {
  input {
    String sample_id
    String sample_type
    String expected_taxon = ""
    File final_calls
    String docker_image = "phemarajata614/afi-terra:0.1.0"
  }

  command <<<
  python3 - <<'PY'
import pandas as pd

calls = pd.read_csv("~{final_calls}", sep="\t")
sample_id = "~{sample_id}"
sample_type = "~{sample_type}"
expected = "~{expected_taxon}".strip()

validation_applicable_types = {"PC_MIX8", "MIXED4", "PC_SINGLE", "CLINICAL", "PC"}

detected = sorted(
    set(
        calls.loc[calls["call"].isin(["Confirmed", "Probable"]), "genus"].astype(str)
    )
)

if sample_type.upper() in validation_applicable_types and expected:
    validation_result = "Concordant" if expected in detected else "Discordant"
else:
    validation_result = "Not_applicable"

out = pd.DataFrame([
    {
        "sample_id": sample_id,
        "sample_type": sample_type,
        "expected_taxon": expected,
        "detected_taxa": ",".join(detected),
        "validation_result": validation_result,
        "validation_control_class": sample_type if sample_type.upper() in {"PC_MIX8", "MIXED4", "PC_SINGLE", "PC"} else "none",
    }
])

out.to_csv("validation_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File validation_summary = "validation_summary.tsv"
  }

  runtime {
    docker: docker_image
  }
}

task SummarizeRoutineTaxa {
  input {
    String sample_id
    String sample_type
    File final_calls
    String docker_image = "phemarajata614/afi-terra:0.1.0"
  }

  command <<<
  python3 - <<'PY'
import pandas as pd

calls = pd.read_csv("~{final_calls}", sep="\t")
sample_id = "~{sample_id}"
sample_type = "~{sample_type}"

detected = sorted(
    set(
        calls.loc[calls["call"].isin(["Confirmed", "Probable"]), "genus"].astype(str)
    )
)

routine_positive_control_types = {"PC_MIX8", "PC"}

out = pd.DataFrame([
    {
        "sample_id": sample_id,
        "sample_type": sample_type,
        "taxa_present": ",".join(detected),
        "n_taxa_present": len(detected),
    "routine_positive_control": "true" if sample_type.upper() in routine_positive_control_types else "false",
    }
])

out.to_csv("routine_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File routine_summary = "routine_summary.tsv"
  }

  runtime {
    docker: docker_image
  }
}
