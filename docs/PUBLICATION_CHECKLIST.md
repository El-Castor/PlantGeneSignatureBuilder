# ✅ PUBLICATION-READY CHECKLIST
## Construction de la liste PCD pour WOT (Waddington-OT)

À inclure dans l'article / méthodes / supplément

---

## 📊 Données & Annotations

- [x] GO annotations issues de eggNOG / HuaCM
- [x] IDs natifs de *Brachypodium distachyon* (Bd21-3)
- [x] **Aucun BLAST inter-espèces utilisé**

---

## 🌱 Définition Biologique

- [x] **Apoptosis remplacée par Programmed Cell Death (PCD)**
- [x] Inclusion explicite des catégories :
  - cell death
  - regulation of cell death
  - oxidative stress / ROS
  - abiotic stress
- [x] Justification biologique spécifique aux plantes

---

## 🔬 Méthodologie

- [x] Mapping protéine → gène déterministe (.v1.2)
- [x] GO filtrés sur **Biological Process uniquement**
- [x] Liste finale filtrée sur les gènes présents dans l'objet single-cell
- [x] Liste utilisée uniquement comme **death score WOT**, pas comme annotation cellulaire

---

## 📝 Reproductibilité

- [x] Script versionné (`create_list_from_GAF.py`)
- [x] Notebook documenté
- [x] GO IDs explicitement listés

---

## 🧬 GO Terms Utilisés

### PCD / Cell Death
- `GO:0012501` — programmed cell death
- `GO:0008219` — cell death
- `GO:0043067` — regulation of programmed cell death
- `GO:0010941` — regulation of cell death

### Stress / ROS
- `GO:0006979` — response to oxidative stress
- `GO:0072593` — reactive oxygen species metabolic process
- `GO:0045454` — cell redox homeostasis
- `GO:0033554` — cellular response to stress
- `GO:0006950` — response to stress
- `GO:0009628` — response to abiotic stimulus

---

## 📐 Schéma Conceptuel

```
GO annotations (eggNOG / HuaCM)
          │
          ▼
Protein-level GO terms
(BdiBd21-3.xGxxxxxxx.1.p)
          │
          ▼
Deterministic ID mapping
(BdiBd21-3.xGxxxxxxx.v1.2)
          │
          ▼
GO filtering (BP only)
PCD + Stress + ROS
          │
          ▼
PCD / death-like gene set
          │
          ▼
WOT death score
(regulates survival / growth term)
```

---

## 🧠 Message Clé

Le score "death" dans WOT capture un **gradient de stress et de mort programmée**,
pas une apoptose animale, et est construit entièrement à partir d'**annotations GO plantes**.

---

## ✅ Protection Reviewer

Cette checklist garantit :
1. **Transparence méthodologique** : GO terms explicites, pas de projection BLAST
2. **Pertinence biologique** : PCD/stress adapté aux plantes
3. **Reproductibilité** : script + notebook + IDs natifs
4. **Utilisation appropriée** : WOT death score, pas annotation cellulaire directe
