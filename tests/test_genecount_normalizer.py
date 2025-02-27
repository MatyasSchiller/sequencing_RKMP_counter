import unittest
import pandas as pd
import numpy as np
import subprocess
import os

class TestGeneExpressionNormalization(unittest.TestCase):

    def setUp(self):
        #Tesztadatok létrehozása.
        self.test_file = "test_input.txt"
        self.output_dir = "./test_output/"
        os.makedirs(self.output_dir, exist_ok=True)

        # Tesztadatok generálása
        test_data = """GeneID\tLength\tSample1\tSample2
gene1\t1000\t50\t100
gene2\t2000\t200\t400
gene3\t1500\t300\t600
"""
        with open(self.test_file, "w") as f:
            f.write(test_data)

    def tearDown(self):
        """Teszt után takarítás."""
        os.remove(self.test_file)
        for file in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, file))
        os.rmdir(self.output_dir)

    def test_rpkm_tpm_calculation(self):
        """Teszteli, hogy az RPKM és TPM számítások működnek-e."""
        df = pd.read_csv(self.test_file, sep="\t")

        gene_lengths = df['Length']
        counts = df.iloc[:, 2:]

        total_mapped_reads = counts.sum(axis=0)
        rpkm = (counts * 1e9) / (gene_lengths.to_numpy()[:, None] * total_mapped_reads.to_numpy()[None, :])
        tpm = (rpkm / rpkm.sum(axis=0)) * 1e6

        self.assertEqual(rpkm.shape, counts.shape)
        self.assertEqual(tpm.shape, counts.shape)
        self.assertTrue(np.all(rpkm >= 0))
        self.assertTrue(np.all(tpm >= 0))

    def test_script_runs(self):
        """Teszteli, hogy a script sikeresen lefut-e."""
        result = subprocess.run(
            ["python3", "genecount_normalizer_main.py", "--input", self.test_file, "--monogram", "AB"],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0)

    def test_missing_input_fails(self):
        """Teszteli, hogy a script hibaüzenetet ad-e, ha nincs input fájl."""
        result = subprocess.run(
            ["python3", "genecount_normalizer_main.py"],
            capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("HIBA: Nem adtál meg input fájlt és nincs stdin bemenet sem!", result.stderr)

if __name__ == "__main__":
    unittest.main()
