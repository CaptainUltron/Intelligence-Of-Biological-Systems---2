"""
AllerScan Backend — app.py
Connects trained_model.pkl to the allergen_ui.html frontend.

Usage:
    pip install flask flask-cors scikit-learn xgboost
    python app.py

Then open http://localhost:5000 in your browser.
"""

import pickle
import math
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

# ── Load model ──────────────────────────────────────────────────
with open("trained_model.pkl", "rb") as f:
    PKL = pickle.load(f)

STAGE_A      = PKL["stage_a_model"]       # XGBClassifier  (binary: allergen vs not)
STAGE_B      = PKL["stage_b_model"]       # RandomForest   (6-class allergen type)
LE6          = PKL["label_encoder_6"]     # decodes Stage B predictions
FEATURE_NAMES = PKL["feature_names"]      # list of 64 feature names (exact order expected by models)

app = Flask(__name__, static_folder=".")
CORS(app)  # allow the HTML file to call this API from any origin


# ── Feature extraction ───────────────────────────────────────────
HYDROPHOBIC_AA = set('AVILMFYW')
CHARGED_AA     = set('DEKRH')
POLAR_AA       = set('STNQ')
AROMATIC_AA    = set('FYW')
ALL_AA         = list('ACDEFGHIKLMNPQRSTVWY')


def extract_features(seq: str) -> dict:
    """
    Extract the 64 features that the model was trained on.
    Physicochemical features + per-amino-acid frequencies + dipeptide frequencies.
    Returns a dict keyed by feature name.
    """
    seq = seq.upper()
    seq = ''.join(c for c in seq if c in set(ALL_AA))
    n = len(seq)
    if n == 0:
        raise ValueError("Empty sequence after cleaning.")

    # --- Per-amino-acid frequencies ---
    aac = {a: seq.count(a) / n for a in ALL_AA}

    # --- Physicochemical group fractions ---
    feats = {
        "hydrophobic":   sum(aac[a] for a in HYDROPHOBIC_AA),
        "charged":       sum(aac[a] for a in CHARGED_AA),
        "polar":         sum(aac[a] for a in POLAR_AA),
        "aromatic":      sum(aac[a] for a in AROMATIC_AA),
        "proline_frac":  aac["P"],
        "cysteine_frac": aac["C"],
        "glycine_frac":  aac["G"],
        "length":        n,
    }

    # --- Per-AA frequencies (aa_X keys) ---
    for a in ALL_AA:
        feats[f"aa_{a}"] = aac[a]

    # --- Dipeptide frequencies (dp_XY keys) ---
    # Only the dipeptides that appear in FEATURE_NAMES are needed.
    dp_features = [f for f in FEATURE_NAMES if f.startswith("dp_")]
    # Build all dipeptide counts first (only once, efficiently)
    dp_counts: dict[str, int] = {}
    for i in range(n - 1):
        dp = seq[i:i+2]
        dp_counts[dp] = dp_counts.get(dp, 0) + 1
    total_dp = max(n - 1, 1)
    for dp_key in dp_features:
        pair = dp_key[3:]   # strip "dp_"
        feats[dp_key] = dp_counts.get(pair, 0) / total_dp

    return feats


def feats_to_vector(feats: dict) -> list:
    """Return feature values in the exact column order the model expects."""
    return [[feats[f] for f in FEATURE_NAMES]]


# ── API endpoint ─────────────────────────────────────────────────
@app.route("/predict", methods=["POST"])
def predict():
    """
    POST /predict
    Body (JSON): { "sequence": "<amino acid string>" }

    Response (JSON):
    {
      "is_allergen": true | false,
      "stage_a_probability": 0.87,
      "allergen_type": "Pollen Allergen" | null,
      "stage_b_probabilities": { "Food Allergen": 0.05, ... } | null
    }
    """
    data = request.get_json(force=True)
    seq  = (data.get("sequence") or "").strip()
    if not seq:
        return jsonify({"error": "No sequence provided."}), 400

    try:
        feats = extract_features(seq)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    X = feats_to_vector(feats)

    # Stage A — binary classification
    stage_a_prob = float(STAGE_A.predict_proba(X)[0][1])   # probability of class 1 (allergen)
    is_allergen  = bool(STAGE_A.predict(X)[0] == 1)

    allergen_type    = None
    stage_b_probs    = None

    # Stage B — allergen type (only if Stage A says allergen)
    if is_allergen:
        stage_b_pred       = STAGE_B.predict(X)[0]
        allergen_type      = LE6.inverse_transform([stage_b_pred])[0]
        stage_b_proba      = STAGE_B.predict_proba(X)[0]
        stage_b_probs      = {
            cls: float(p)
            for cls, p in zip(LE6.classes_, stage_b_proba)
        }

    return jsonify({
        "is_allergen":          is_allergen,
        "stage_a_probability":  stage_a_prob,
        "allergen_type":        allergen_type,
        "stage_b_probabilities": stage_b_probs,
    })


# ── Serve the HTML UI ────────────────────────────────────────────
@app.route("/")
def index():
    return send_from_directory(".", "allergen_ui.html")


if __name__ == "__main__":
    print("AllerScan backend running — open http://localhost:5000")
    app.run(debug=True, port=5000)
