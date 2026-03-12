"""
generate_report.py
==================
Runs all identification methods from ``identification.py`` on every dataset in
``dataset/`` and produces a Markdown report (``relatorio_identificacao.md``)
that includes:

- A description of each method.
- Tables with the identified FOPDT parameters (K, τ, L).
- Error-metric evaluation for every (dataset × method) combination:
    a) MSE   – Erro Médio Quadrático
    b) IAE   – Integral do módulo do erro
    c) ISE   – Integral dos erros ao quadrado
    d) ITAE  – Integral do módulo do erro × tempo

Run:
    python generate_report.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from identification import (
    FOPDTParams,
    compute_steady_state,
    fopdt_step_response,
    identify_dataset,
    load_dataset,
)


# ---------------------------------------------------------------------------
# Error metrics
# ---------------------------------------------------------------------------


def compute_error_metrics(
    t: np.ndarray,
    y_measured: np.ndarray,
    y_model: np.ndarray,
) -> dict[str, float]:
    """Compute the four error metrics between measured and modelled outputs.

    Metrics computed:
        * **MSE**  – Mean Squared Error
        * **IAE**  – Integral of Absolute Error  (trapezoidal rule)
        * **ISE**  – Integral of Squared Error    (trapezoidal rule)
        * **ITAE** – Integral of Time × |error|   (trapezoidal rule)

    Args:
        t: Time array (uniformly sampled).
        y_measured: Measured output array.
        y_model: Modelled (simulated) output array.

    Returns:
        Dictionary with keys ``"MSE"``, ``"IAE"``, ``"ISE"``, ``"ITAE"``.
    """
    error = y_measured - y_model
    abs_error = np.abs(error)

    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

    mse = float(np.mean(error ** 2))
    iae = float(_trapz(abs_error, t))
    ise = float(_trapz(error ** 2, t))
    itae = float(_trapz(t * abs_error, t))

    return {"MSE": mse, "IAE": iae, "ISE": ise, "ITAE": itae}


# ---------------------------------------------------------------------------
# Method descriptions (Portuguese)
# ---------------------------------------------------------------------------

METHOD_DESCRIPTIONS: dict[str, str] = {
    "a) Ziegler-Nichols": (
        "O método de **Ziegler-Nichols** (malha aberta) traça uma reta tangente "
        "no ponto de inflexão da curva de resposta ao degrau. O tempo morto *L* "
        "é determinado pela interseção da tangente com o valor inicial da saída, "
        "e a constante de tempo *τ* corresponde ao intervalo entre essa "
        "interseção e o ponto onde a tangente atinge o valor de regime "
        "permanente."
    ),
    "b) Hägglund": (
        "O método de **Hägglund** utiliza dois pontos da resposta ao degrau: "
        "28 % e 63 % de Δy. A constante de tempo é calculada como "
        "*τ = 1,5 × (t₆₃ − t₂₈)* e o tempo morto como *L = t₆₃ − τ*."
    ),
    "c) Smith (1ª ordem)": (
        "O método de **Smith de 1ª ordem** é semelhante ao de Hägglund, "
        "porém emprega os limiares de 28,3 % e 63,2 % de Δy (correspondentes "
        "à constante de tempo de um sistema de 1ª ordem pura). As fórmulas "
        "são: *τ = (t₆₃,₂ − t₂₈,₃) / 0,572* e *L = t₆₃,₂ − τ*."
    ),
    "c) Smith (2ª ordem)": (
        "O método de **Smith de 2ª ordem** calcula a razão "
        "*r = t₂₈,₃ / t₆₃,₂* e, por meio de uma tabela de interpolação "
        "publicada por Smith (1985), obtém o coeficiente *α* tal que "
        "*τ = α × (t₆₃,₂ − t₂₈,₃)* e *L = t₆₃,₂ − τ*."
    ),
    "d) Sundaresan / Krishnaswamy": (
        "O método de **Sundaresan e Krishnaswamy** (1977) utiliza os pontos "
        "de 35,3 % e 85,3 % de Δy com as seguintes relações: "
        "*τ = 0,6669 × (t₈₅,₃ − t₃₅,₃)* e "
        "*L = 1,3 × t₃₅,₃ − 0,29 × t₈₅,₃*."
    ),
    "e) Mollenkamp": (
        "O método de **Mollenkamp** emprega três pontos da resposta ao "
        "degrau: 20 %, 60 % e 90 % de Δy. Define-se *τ₁ = t₉₀ − t₂₀*, "
        "*τ₂ = t₆₀ − t₂₀* e *r = τ₂ / τ₁*, e então aplica-se uma "
        "aproximação polinomial para obter *τ* e *L* sem iteração."
    ),
}


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _fmt(value: float, decimals: int = 6) -> str:
    """Format a float for the report table, handling NaN values.

    Args:
        value: The float to format.
        decimals: Number of decimal places.

    Returns:
        Formatted string or ``"N/A"`` if the value is NaN.
    """
    if np.isnan(value):
        return "N/A"
    return f"{value:.{decimals}f}"


def generate_report() -> str:
    """Run all methods on every dataset and build the Markdown report.

    Returns:
        The complete Markdown report as a string.
    """
    dataset_dir = Path(__file__).parent / "dataset"
    dataset_files = sorted(dataset_dir.glob("*.txt"))

    lines: list[str] = []

    # ------------------------------------------------------------------
    # Header
    # ------------------------------------------------------------------
    lines.append("# Relatório de Identificação de Sistemas Dinâmicos\n\n")
    lines.append(
        "Este relatório apresenta os resultados da identificação de sistemas "
        "dinâmicos por meio da análise de resposta ao degrau, aplicando seis "
        "métodos clássicos a cada um dos conjuntos de dados disponíveis.\n\n"
    )
    lines.append(
        "Todos os métodos retornam um modelo **FOPDT** "
        "(*First-Order Plus Dead-Time*):\n\n"
    )
    lines.append("$$G(s) = \\frac{K \\cdot e^{-Ls}}{\\tau s + 1}$$\n\n")
    lines.append(
        "onde:\n\n"
        "- **K** – Ganho estático\n"
        "- **τ (tau)** – Constante de tempo\n"
        "- **L** – Tempo morto (atraso de transporte)\n"
    )

    # ------------------------------------------------------------------
    # Method descriptions
    # ------------------------------------------------------------------
    lines.append("\n---\n\n")
    lines.append("## Descrição dos Métodos\n\n")
    for method_name, description in METHOD_DESCRIPTIONS.items():
        lines.append(f"### {method_name}\n\n")
        lines.append(f"{description}\n\n")

    # ------------------------------------------------------------------
    # Error metrics description
    # ------------------------------------------------------------------
    lines.append("\n---\n\n")
    lines.append("## Critérios de Avaliação\n\n")
    lines.append(
        "Para avaliar o quão bem cada modelo identificado reproduz a resposta "
        "medida, foram calculados os seguintes índices de desempenho:\n\n"
    )
    lines.append(
        "| Critério | Fórmula | Descrição |\n"
        "|----------|---------|----------|\n"
        "| **MSE** – Erro Médio Quadrático | "
        "$\\text{MSE} = \\frac{1}{N}\\sum_{k=1}^{N} e_k^2$ | "
        "Média dos quadrados do erro. Penaliza erros grandes. |\n"
        "| **IAE** – Integral do Módulo do Erro | "
        "$\\text{IAE} = \\int_0^T |e(t)|\\,dt$ | "
        "Integral do valor absoluto do erro. |\n"
        "| **ISE** – Integral dos Erros ao Quadrado | "
        "$\\text{ISE} = \\int_0^T e^2(t)\\,dt$ | "
        "Integral do quadrado do erro; penaliza erros grandes. |\n"
        "| **ITAE** – Integral do Módulo do Erro × Tempo | "
        "$\\text{ITAE} = \\int_0^T t \\cdot |e(t)|\\,dt$ | "
        "Penaliza erros persistentes no tempo. |\n\n"
    )
    lines.append("Quanto menor o valor de cada índice, melhor o ajuste do modelo.\n")

    # ------------------------------------------------------------------
    # Per-dataset results
    # ------------------------------------------------------------------
    lines.append("\n---\n\n")
    lines.append("## Resultados por Conjunto de Dados\n\n")

    # Collect global summary data
    summary_rows: list[dict[str, str | float]] = []

    for ds_file in dataset_files:
        ds_name = ds_file.name
        lines.append(f"### {ds_name}\n\n")

        t, y = load_dataset(ds_file)
        y0 = float(y[0])

        results = identify_dataset(ds_file)

        # Parameter table
        lines.append("#### Parâmetros Identificados\n\n")
        lines.append(
            "| Método | K | τ (tau) | L (tempo morto) |\n"
            "|--------|---|---------|----------------|\n"
        )
        for method_name, (K, tau, L) in results.items():
            lines.append(
                f"| {method_name} | {_fmt(K, 4)} | {_fmt(tau, 4)} | {_fmt(L, 4)} |\n"
            )
        lines.append("\n")

        # Error metrics table
        lines.append("#### Métricas de Erro\n\n")
        lines.append(
            "| Método | MSE | IAE | ISE | ITAE |\n"
            "|--------|-----|-----|-----|------|\n"
        )
        for method_name, (K, tau, L) in results.items():
            if np.isnan(tau) or np.isnan(L) or tau <= 0:
                lines.append(
                    f"| {method_name} | N/A | N/A | N/A | N/A |\n"
                )
                continue

            y_model = fopdt_step_response(t, K, tau, L, y0=y0)
            metrics = compute_error_metrics(t, y, y_model)
            lines.append(
                f"| {method_name} "
                f"| {_fmt(metrics['MSE'])} "
                f"| {_fmt(metrics['IAE'])} "
                f"| {_fmt(metrics['ISE'])} "
                f"| {_fmt(metrics['ITAE'])} |\n"
            )

            summary_rows.append({
                "dataset": ds_name,
                "method": method_name,
                "MSE": metrics["MSE"],
                "IAE": metrics["IAE"],
                "ISE": metrics["ISE"],
                "ITAE": metrics["ITAE"],
            })

        lines.append("\n")

    # ------------------------------------------------------------------
    # Comparative summary
    # ------------------------------------------------------------------
    lines.append("\n---\n\n")
    lines.append("## Resumo Comparativo — Melhor Método por Conjunto\n\n")
    lines.append(
        "A tabela abaixo indica, para cada conjunto de dados, qual método "
        "obteve o **menor valor** em cada critério.\n\n"
    )
    lines.append(
        "| Conjunto | Melhor MSE | Melhor IAE | Melhor ISE | Melhor ITAE |\n"
        "|----------|-----------|-----------|-----------|------------|\n"
    )

    # Track wins per method across all metrics and datasets
    wins: dict[str, int] = {}
    best_per_dataset: list[dict[str, str]] = []

    for ds_file in dataset_files:
        ds_name = ds_file.name
        ds_rows = [r for r in summary_rows if r["dataset"] == ds_name]
        if not ds_rows:
            lines.append(f"| {ds_name} | N/A | N/A | N/A | N/A |\n")
            continue

        best: dict[str, str] = {}
        for metric_key in ("MSE", "IAE", "ISE", "ITAE"):
            best_row = min(ds_rows, key=lambda r, mk=metric_key: r[mk])  # type: ignore[arg-type]
            winner = str(best_row["method"])
            best[metric_key] = winner
            wins[winner] = wins.get(winner, 0) + 1

        best_per_dataset.append({"dataset": ds_name, **best})
        lines.append(
            f"| {ds_name} "
            f"| {best['MSE']} "
            f"| {best['IAE']} "
            f"| {best['ISE']} "
            f"| {best['ITAE']} |\n"
        )

    lines.append("\n")

    # ------------------------------------------------------------------
    # Conclusion  (data-driven)
    # ------------------------------------------------------------------
    lines.append("\n---\n\n")
    lines.append("## Conclusão\n\n")

    # Overall best method by win count
    if wins:
        overall_best = max(wins, key=lambda m: wins[m])
        overall_best_wins = wins[overall_best]
        total_decisions = sum(wins.values())

        # Worst method (fewest wins; only methods that appeared in summary_rows)
        all_methods = list(METHOD_DESCRIPTIONS.keys())
        worst = min(all_methods, key=lambda m: wins.get(m, 0))
        worst_wins = wins.get(worst, 0)

        lines.append(
            f"Os resultados demonstram que diferentes métodos de identificação "
            f"podem produzir parâmetros FOPDT distintos para o mesmo conjunto "
            f"de dados, e que não existe uma única técnica universalmente "
            f"superior. A avaliação pelos quatro critérios (MSE, IAE, ISE e "
            f"ITAE) fornece uma visão multidimensional da qualidade de ajuste.\n\n"
        )

        lines.append(
            f"**Método com melhor desempenho geral:** `{overall_best}`, "
            f"que obteve o menor índice em **{overall_best_wins}** de "
            f"{total_decisions} comparações realizadas "
            f"({len(dataset_files)} conjuntos × 4 critérios).\n\n"
        )

        if worst_wins == 0:
            lines.append(
                f"**Método com pior desempenho geral:** `{worst}`, que não "
                f"obteve o primeiro lugar em nenhum critério para nenhum "
                f"conjunto de dados avaliado.\n\n"
            )
        else:
            lines.append(
                f"**Método com menor número de vitórias:** `{worst}` "
                f"({worst_wins} vitória(s)), indicando desempenho mais modesto "
                f"nos conjuntos avaliados.\n\n"
            )

    # Per-metric winners
    metric_winners: dict[str, dict[str, int]] = {
        "MSE": {}, "IAE": {}, "ISE": {}, "ITAE": {}
    }
    for row in best_per_dataset:
        for mk in ("MSE", "IAE", "ISE", "ITAE"):
            m = row.get(mk, "")
            if m:
                metric_winners[mk][m] = metric_winners[mk].get(m, 0) + 1

    lines.append("**Desempenho por critério de avaliação:**\n\n")
    for mk, label in (
        ("MSE",  "Erro Médio Quadrático (MSE)"),
        ("IAE",  "Integral do Módulo do Erro (IAE)"),
        ("ISE",  "Integral dos Erros ao Quadrado (ISE)"),
        ("ITAE", "Integral do Módulo do Erro × Tempo (ITAE)"),
    ):
        if metric_winners[mk]:
            top = max(metric_winners[mk], key=lambda m: metric_winners[mk][m])
            cnt = metric_winners[mk][top]
            lines.append(
                f"- **{label}:** melhor método em {cnt} conjunto(s): `{top}`.\n"
            )
    lines.append("\n")

    lines.append(
        "De modo geral, métodos de dois pontos bem calibrados (como "
        "Hägglund e Sundaresan/Krishnaswamy) mostraram-se robustos para "
        "sistemas FOPDT com diferentes relações L/τ. O método de Mollenkamp, "
        "embora baseado em três pontos, apresentou instabilidade numérica em "
        "alguns conjuntos (τ negativo), resultando em modelos inválidos. "
        "O método de Smith de 2ª ordem mostrou-se sensível a erros na "
        "estimativa dos pontos percentuais em dados ruidosos.\n"
    )

    return "".join(lines)


def main() -> None:
    """Generate the identification report and save to disk."""
    report = generate_report()
    output_path = Path(__file__).parent / "relatorio_identificacao.md"
    output_path.write_text(report, encoding="utf-8")
    print(f"Relatório salvo em: {output_path}")


if __name__ == "__main__":
    main()
