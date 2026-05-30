import os
import json
import subprocess
import yaml

# =====================================================================
# 1. SETUP PARAMETERS & PATHS
# =====================================================================
# Isolate Domain III of the EGFR extracellular domain (PDB 8HGO Chain A)
EGFR_DOMAIN_III_SEQ = "LEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFM"

# Path to the JSON containing your variants pool and the template PDB
JSON_INPUT_FILE = "./data_inputs/random_variants.json"
TEMPLATE_PDB_FILE = "./data_inputs/1YY9.pdb"

# Working directory to save files
OUTPUT_ROOT = "./boltz_results"
YAML_DIR = os.path.join(OUTPUT_ROOT, "input_yaml_configs")

os.makedirs(YAML_DIR, exist_ok=True)

# Helper custom dumper to handle block styles natively
class ExplicitDumper(yaml.SafeDumper):
    pass

# Ensure lists are formatted cleanly inside brackets [A, B] as requested
def list_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
ExplicitDumper.add_representer(list, list_representer)

# =====================================================================
# 2. READ VARIANT POOL FROM JSON
# =====================================================================
if not os.path.exists(JSON_INPUT_FILE):
    raise FileNotFoundError(f"Missing required input file: {JSON_INPUT_FILE}. Ensure you generated your variants pool first.")

with open(JSON_INPUT_FILE, "r") as f:
    variant_pool = json.load(f)

print(f"📋 Loaded {len(variant_pool)} variants for processing.")

# =====================================================================
# 3. GENERATE THE CORRECT BOLTZ FORMAT YAML & SUBMIT FOR CO-FOLDING
# =====================================================================
for variant in variant_pool:
    v_id = variant["id"]
    v_seq = variant["seq"]
    
    print(f"\n--- Processing Variant: {v_id} ---")
    
    # Constructing the exact data block template requested
    boltz_config = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": v_seq,
                    "msa": "empty"
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": EGFR_DOMAIN_III_SEQ,
                    "msa": "empty"
                }
            }
        ],
        "templates": [
            {
                "pdb": TEMPLATE_PDB_FILE,
                "chain_id": ["A", "B"],
                "force": True,
                "threshold": 3.0
            }
        ]
    }
    
    # Write the YAML configuration out
    yaml_filename = f"{v_id}_config.yaml"
    yaml_path = os.path.join(YAML_DIR, yaml_filename)
    
    with open(yaml_path, 'w') as file:
        yaml.dump(boltz_config, file, Dumper=ExplicitDumper, default_flow_style=False, sort_keys=False)
        
    print(f"💾 Generated YAML configuration schema at: {yaml_path}")
    
    # Establish sub-directories for each variant result output path
    variant_out_dir = os.path.join(OUTPUT_ROOT, f"{v_id}_run")
    
    # Build CLI command structure
    # --use_msa_server is omitted here because 'msa: empty' is explicitly requested
    cmd = [
        "boltz", "predict", yaml_path,
        "--out_dir", variant_out_dir,
        "--recycling_steps", "1",  # CRITICAL: Drops recycling iterations from 3 to 1 for Mac optimization
        "--num_workers", "1",
        '--use_msa_server'       
    ]

    print(f"🛸 Submitting task execution: `{' '.join(cmd)}`")

    subprocess.run(cmd, check=True)
