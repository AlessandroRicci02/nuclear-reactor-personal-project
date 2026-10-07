"""Write a compact human-readable snapshot from completed search assessments."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def fmt_k(result):
    if result is None:
        return "—"
    return f"{result['k_eff']:.6f} ± {result['sigma_k']:.6f}"


def fmt_alpha(result):
    if result is None:
        return "—"
    return (f"{result['alpha_pcm_K']:+.3f} ± {result['sigma_pcm_K']:.3f} "
            f"(limite superiore 95%: {result['upper_95_pcm_K']:+.3f})")


def main():
    assessments = []
    for path in (ROOT / "campaign").glob("*/assessment.json"):
        assessments.append(json.loads(path.read_text()))
    assessments.sort(key=lambda item: item["index"])
    status_path = ROOT / "status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    lines = ["RICERCA NEUTRONICA YH2 PER NEP — STATO CORRENTE", "",
             "Criteri: k_OUT>1, k_IN<1 e coefficienti di combustibile, "
             "moderatore e totale <0, con margini statistici indicati.",
             "I coefficienti di moderatore e totale includono la stima "
             "conservativa di dilatazione YH2 beta=2.4e-5/K.",
             "Le incertezze riportate sono Monte Carlo (1 sigma): non includono "
             "incertezza delle sezioni d'urto, della legge S(a,b) o della "
             "dilatazione di componenti reali.",
             "Il profilo He nominale 800→1150 K è imposto, non verificato "
             "termoidraulicamente.", "",
             f"Candidati completati: {len(assessments)}",
             f"Processo: {status.get('status', 'ignoto')}; "
             f"candidato attuale: {status.get('candidate', '—')}; "
             f"fase: {status.get('stage', '—')}; caso: {status.get('case', '—')}", ""]
    for item in assessments:
        final = item.get("confirmation") or item.get("pilot") or {}
        beta_max = item.get("beta_max") or {}
        results = final.get("results", {})
        coefficients = final.get("coefficients", {})
        lines.extend([
            f"[{item['index']:02d}] {item['candidate']} — {item['status']}",
            f"  motivo: {beta_max.get('reason', final.get('reason', item.get('reason', '—')))}",
            f"  k_OUT: {fmt_k(results.get('out'))}; "
            f"k_IN: {fmt_k(results.get('in'))}",
            f"  alpha_moderatore [pcm/K]: {fmt_alpha(coefficients.get('moderator'))}",
            f"  alpha_combustibile [pcm/K]: {fmt_alpha(coefficients.get('fuel'))}",
            f"  alpha_totale [pcm/K]: {fmt_alpha(coefficients.get('total'))}",
        ])
        if beta_max:
            end_coeff = beta_max.get("coefficients", {})
            lines.extend([
                f"  alpha_moderatore beta_max [pcm/K]: "
                f"{fmt_alpha(end_coeff.get('moderator'))}",
                f"  alpha_totale beta_max [pcm/K]: "
                f"{fmt_alpha(end_coeff.get('total'))}",
            ])
        lines.append("")
    (ROOT / "progress.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
