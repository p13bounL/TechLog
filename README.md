scRNA-seq Pipeline App (Streamlit)

Διαδραστική εφαρμογή Streamlit για εκτέλεση βασικού scRNA-seq pipeline:
- Μετατροπή 10x mtx → `.h5ad`
- Προεπεξεργασία / QC
- Συνένωση δειγμάτων
- Ενσωμάτωση (Scanorama)
- Annotation (Leiden-based)
- DEG Analysis (2 ομάδες)
- Volcano & Expression plots
- Export ρυθμίσεων/plots/DEG

Το project είναι **Dockerized**, τρέχει σε Windows/macOS/Linux μόνο με Docker.

---

## 🔧 Απαιτήσεις

- [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- 8–16 GB RAM (ανάλογα το μέγεθος δεδομένων)

---

## 🚀 Γρήγορη εκκίνηση

Κλωνοποίησε ή κατέβασε το repo, άνοιξε τερματικό **στον φάκελο του repo** και:

```bash
# Build image
docker build -t scapp:latest .

# Τρέξε app (port 8501) και κάνε mount τα δεδομένα σου
# Linux/macOS:
docker run --rm -p 8501:8501 -v "$PWD/data:/app/data" scapp:latest

# Windows PowerShell:
docker run --rm -p 8501:8501 -v "${PWD}\data:/app/data" scapp:latest