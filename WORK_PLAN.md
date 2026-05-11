WORK PLAN — qewan
=================

Resumo do que já foi feito
--------------------------
- Código atualizado em `qewan/` (vários módulos): CLI e geradores.
- CLI (`qewan/cli.py`):
  - Adicionadas flags: `--write-mmn/--no-write-mmn`, `--write-amn/--no-write-amn`, `--write-unk/--no-write-unk` em `run_pw2wannier`.
  - `--auto-projections/--no-auto-projections` exposto em `run_pw2wannier`.
  - `--include-projections/--no-include-projections` adicionado ao `run_wan`.
  - `--conventional-cell/--no-conventional-cell` (default: False) adicionado a `run_scf`, `run_bands`, `run_nscf`, `run_wan` — controla se usa a célula convencional do ASE ou reduz por simetria.
- IO / geradores (`qewan/io.py`):
  - `atoms_to_pw_input`: agora escreve `ATOMIC_POSITIONS crystal` (coordenadas fracionárias) para SCF/BANDS/NSCF.
  - `generate_pw2wannier_input`: `atom_proj` é escrito apenas quando `auto_projections=True`; aceita flags `write_mmn/write_amn/write_unk`.
  - `generate_wannier_win`: escreve `mp_grid` imediatamente acima de `begin kpoints` caso aplicável; o bloco `begin kpoints` contém apenas coordenadas (sem pesos).
  - `_squeeze_blank_lines`: função utilitária para colapsar quebras de linha em excesso nas saídas.
  - `reduce_atoms_by_symmetry(atoms, ...)`: nova função que usa `spglib` (se instalado) para retornar célula primitiva ou representantes únicos e mapping de equivalência.
- Documentação:
  - `README.md` atualizado: indica que `ATOMIC_POSITIONS` são fracionárias e documenta `--conventional-cell` e requisito/opção `spglib`.
- Limpeza:
  - A pasta `demo_run/` criada nos testes foi removida conforme pedido.
- Testes manuais executados:
  - Geração de SCF/BANDS/NSCF e `.win` para `Gd.cif` e `GdCo2.cif` (manualmente, verificados); verificação de que `atom_proj` aparece/é omitido conforme `auto_projections`.
  - Verificado que redução por simetria transforma 24 átomos → 6 (primitive) no `GdCo2.cif` de exemplo.

Observações técnicas importantes
--------------------------------
- A redução por simetria usa `spglib`. Se `spglib` não estiver instalado, a redução é ignorada — instrução para instalar: `pip install spglib`.
- `ATOMIC_POSITIONS` são sempre escritas como `crystal` (coordenadas fracionárias), independentemente de estar usando célula convencional ou reduzida.

Plano de trabalho sugerido (itens possíveis a implementar)
---------------------------------------------------------
Ordem proposta por prioridade (curto → longo prazo):

1) Expor `pw2wannier` options em `run_wan` (alta prioridade)
   - Tornar `run_wan` capaz de gerar `pw2wannier.in` diretamente com as mesmas flags: `--seedname`, `--auto-projections/--no-auto-projections`, `--write-mmn`, `--write-amn`, `--write-unk`.
   - Motivo: fluxo end-to-end mais simples (do CIF → .win + pw2wannier).
   - Estimativa: 30–60 min.

2) `--reduce-mode` (média prioridade)
   - Adicionar `--reduce-mode=primitive|representatives|none` para controlar comportamento de redução por simetria:
     - `primitive` (atual): retorna célula primitiva via `spglib.find_primitive`.
     - `representatives`: retornar apenas representantes únicos por classe (`equivalent_atoms` → uma posição por classe).
     - `none`: não reduzir.
   - Motivo: casos onde se quer apenas uma lista de representantes (p.ex. construção de projeções) vs célula primitiva completa.
   - Estimativa: 30–60 min.

3) Expor `write_mmn/write_amn/write_unk` e `auto_projections` também em `run_wan` (média prioridade)
   - Motivo: consistência CLI e evitar necessidade de rodar `run_pw2wannier` separadamente.
   - Estimativa: 15–30 min.

4) Testes automatizados (alta prioridade)
   - Adicionar `pytest` com fixtures CIFs pequenas, testar:
     - `atoms_to_pw_input` formato (presença de `ATOMIC_POSITIONS crystal`).
     - `generate_wannier_win` (mp_grid, ausência/presença de `projections`).
     - `generate_pw2wannier_input` (atom_proj condicional, flags write_*).
   - Simular ausência de `spglib` para teste de fallback.
   - Estimativa: 2–4 horas.

5) Quality / CI / formatting (média)
   - Rodar `black`/`ruff` e adicionar GitHub Action para rodar `pytest` + `black --check`.
   - Fazer commit e abrir PR.
   - Estimativa: 1–2 horas.

6) Documentação e exemplos (curto)
   - Adicionar exemplos no `README` para:
     - `--conventional-cell` e `--reduce-mode`.
     - Exemplo end-to-end curto com `run_wan` → `run_pw2wannier` (quando exposto).
   - Estimativa: 15–30 min.

7) Integração com busca automática de pseudopotenciais (médio)
   - Permitir especificar `--pseudos-dir` ou procurar `PSEUDO_DIR` configurado.
   - Estimativa: 30–60 min.

8) Melhorias finas (opcionais)
   - Validadores de input CLI (ex.: formato correto de `--kmesh`).
   - Opção para escrever `ATOMIC_POSITIONS angstrom` quando o usuário explicitamente pedir (atualmente sempre `crystal`).
   - Support for spin/polarization flags in pw.x generation.

Arquivo de status/local
------------------------
- Este arquivo: `WORK_PLAN.md` foi salvo na raiz do repositório para referência futura.
- O TODO interno também foi atualizado via ferramenta de gerenciamento (`manage_todo_list`).

Como continuar daqui (comandos úteis)
-------------------------------------
- Rodar um SCF com célula reduzida (padrão):

```bash
python -m qewan.cli run_scf path/to/structure.cif --outdir run_scf --kpoints "6 6 6"
```

- Rodar um SCF usando a célula convencional do CIF (sem redução):

```bash
python -m qewan.cli run_scf path/to/structure.cif --outdir run_scf --kpoints "6 6 6" --conventional-cell
```

- Gerar `.win` com autoproj (padrão):

```bash
python -m qewan.cli run_wan path/to/structure.cif --outdir run_wan --nscf-dir run_nscf --seedname myseed --auto-projections
```

- Gerar `pw2wannier.in` (agora exposto em `run_pw2wannier`):

```bash
python -m qewan.cli run_pw2wannier path/to/structure.cif --outdir run_pw2wannier --prefix myprefix --seedname myseed --write-mmn --write-amn --no-write-unk --auto-projections
```

Notas finais
------------
- Se quiser, eu implemento o item (1) agora (expor `pw2wannier` options diretamente em `run_wan`) e atualizo os exemplos do `README`, ou eu posso abrir um branch/commit com tudo que foi feito e aplicar formatação (`black`) antes de preparar um PR.
- Diga qual opção prefere e eu executo em seguida.
