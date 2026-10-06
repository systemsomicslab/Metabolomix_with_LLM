import json
import unittest

import server
from metabolomix.core import session_state


class TestVerifyHasIdentityBlock(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_dossier_includes_identity_normalization(self):
        feat = {"id": 1, "name": "PC 34:1", "ontology": "PC",
                "formula": "C42H82NO8P", "adduct": "[M+HCOO]-",
                "m/z": 850.0, "ion_mode": "Negative", "has_msms": True,
                "time": {"rt": 5.0}}
        vocab = {}
        dossier = server._build_verification_dossier(feat, vocab)
        self.assertIn("identity_normalization", dossier)
        self.assertIn("msi", dossier["identity_normalization"])


if __name__ == "__main__":
    unittest.main()
