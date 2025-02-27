import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os
import sys
import argparse
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages

# Argumentumok beolvasása
parser = argparse.ArgumentParser(description="Gene expression normalization (RPKM, TPM) and visualization")
parser.add_argument("--input", type=str, help="Path to input file (default: read from stdin)", default=None)
parser.add_argument("--monogram", type=str, help="User monogram (default: 'SM')", default="SM")
args = parser.parse_args()

# Failsafe: se input fájl, se stdin nincs -> hiba
if not args.input and sys.stdin.isatty():
    print("HIBA: Nem adtál meg input fájlt és nincs stdin bemenet sem!", file=sys.stderr)
    sys.exit(1)

# Input fájl betöltése
if args.input:
    file_path = args.input
    df = pd.read_csv(file_path, sep="\t", comment='#')
else:
    df = pd.read_csv(sys.stdin, sep="\t", comment='#')

monogram = args.monogram
now = datetime.now().strftime("%Y%m%d_%H%M")
output_dir = f"./20221124_skill_survey__rscript/data/{now}/"
os.makedirs(output_dir, exist_ok=True)
log_file_path = f"{output_dir}/skill_survey__abundance_distribution_{now}_1{monogram}.log"

# Korábbi log handlerek eltávolítása
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(
    filename=log_file_path,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.info("Script elindult.")
if args.input:
    logging.info(f"Adatok beolvasása a következő fájlból: {file_path}")
else:
    logging.info("Adatok beolvasása stdin-ből.")

# 2. RPKM számítása
gene_lengths = df['Length']
counts = df.iloc[:, 6:]

total_mapped_reads = counts.sum(axis=0)
logging.info("RPKM számítás elkezdődött...")
rpkm = (counts * 1e9) / (gene_lengths.to_numpy()[:, None] * total_mapped_reads.to_numpy()[None, :])
logging.info("RPKM számítás befejeződött.")

# 3. TPM számítása
logging.info("TPM számítás elkezdődött...")
tpm = (rpkm / rpkm.sum(axis=0)) * 1e6
logging.info("TPM számítás sikeresen befejeződött.")

# 4. CSV fájlok mentése
rpkm.to_csv(f"{output_dir}/skill_survey__rpkm_{now}_1.csv", index=False)
tpm.to_csv(f"{output_dir}/skill_survey__tpm_{now}_1.csv", index=False)

# 5. Log-transzformáció és ábrázolás
output_pdf_path = f"{output_dir}/skill_survey__abundance_distribution_{now}_1{monogram}.pdf"
with PdfPages(output_pdf_path) as pdf:
    plt.figure()
    for col in counts.columns:
        sns.kdeplot(np.log1p(counts[col]), label=col)
    plt.title("Density of raw gene counts")
    plt.xlabel("Log-transformed gene count")
    plt.ylabel("Density")
    plt.legend()
    pdf.savefig()
    plt.close()

    plt.figure()
    for col in tpm.columns:
        sns.kdeplot(np.log1p(tpm[col]), label=col)
    plt.title("Density of TPM-normalized gene counts")
    plt.xlabel("Log-transformed gene count")
    plt.ylabel("Density")
    plt.legend()
    pdf.savefig()
    plt.close()

logging.info("Sűrűségfüggvények ábrázolása befejeződött.")
logging.info(f"Grafikonok elmentve ide: {output_pdf_path}")
logging.info("Script sikeresen lefutott.")
