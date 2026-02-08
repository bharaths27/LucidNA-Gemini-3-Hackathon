import sys
import subprocess
import pkg_resources

# --- AUTO-FIX: Install missing libraries at runtime ---
required = {'py3dmol', 'stmol', 'google-generativeai'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print(f"Installing missing libraries: {missing}")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
    print("Libraries installed! Please restart the app if it doesn't load automatically.")
# -------------------------------------------------------


import streamlit as st
import requests
import time
import google.generativeai as genai
from stmol import showmol
import py3Dmol


# --- CONFIGURATION & SECRETS ---
st.set_page_config(page_title="LucidDNA Sentry", layout="wide", page_icon="🧬")

# Try to load API key from secrets, or use a placeholder if missing
# Make sure to create .streamlit/secrets.toml with: GEMINI_API_KEY = "your_key"
api_key = st.secrets.get("GEMINI_API_KEY", None)
if api_key:
    genai.configure(api_key=api_key)

# --- MOCK DATA (The "Monke" Strategy) ---
DEMO_DATA = {
    "BRCA1": {
        "variant": "c.5266dupC (p.Gln1756Profs)",
        "risk": "High (Pathogenic)",
        "pdb_id": "1JM7",  # Real PDB ID for BRCA1 BRCT domain
        "mutation_resi": "50", # Relative position in this PDB fragment
        "description": "Frameshift mutation in the BRCT domain, disrupting DNA repair binding.",
        "tissue": "Breast / Ovarian",
        "score": -14.2
    },
    "TP53": {
        "variant": "c.524G>A (p.Arg175His)",
        "risk": "Critical (Pathogenic)",
        "pdb_id": "1TSR",  # Real PDB ID for TP53 DNA binding domain
        "mutation_resi": "175",
        "description": "Hotspot mutation in the DNA-binding domain. Destabilizes the protein core.",
        "tissue": "Pan-Tissue",
        "score": -18.5
    },
    "CFTR": {
        "variant": "c.1521_1523delCTT (p.Phe508del)",
        "risk": "High (Pathogenic)",
        "pdb_id": "1XMI", 
        "mutation_resi": "508",
        "description": "Deletion of Phenylalanine at 508 causing protein misfolding.",
        "tissue": "Lung / Pancreas",
        "score": -10.1
    }
}

# --- HELPER FUNCTIONS ---
def render_protein(pdb_id, resi_to_highlight):
    """Fetches PDB data server-side and renders it to avoid CORS issues."""
    try:
        # Fetch directly from RCSB PDB
        url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code != 200:
            return None, f"Error fetching PDB {pdb_id}"
        
        pdb_data = response.text
        
        # Configure the 3D View
        view = py3Dmol.view(width=500, height=400)
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Highlight the mutation site (Red Sphere)
        view.addStyle({'resi': resi_to_highlight}, {'sphere': {'color': 'red', 'opacity': 0.8}})
        view.addLabel(f"Mutation: {resi_to_highlight}", 
                     {'position': {'resi': resi_to_highlight}, 
                      'backgroundColor': 'black', 
                      'fontColor': 'white'})
        
        view.zoomTo()
        return view, None
    except Exception as e:
        return None, str(e)

def get_gemini_interpretation(gene, data):
    """Generates a clinical explanation using Gemini."""
    if not api_key:
        return "Gemini API Key missing. Please check .streamlit/secrets.toml."
    
    try:
        model = genai.GenerativeModel('gemini-2.5-flash')
        prompt = f"""
        You are an expert geneticist. Interpret this variant for a patient report.
        
        Gene: {gene}
        Variant: {data['variant']}
        Risk: {data['risk']}
        Structural Impact: {data['description']}
        
        Explain WHY this specific mutation is dangerous in plain English. 
        Reference the protein structure (e.g., "The mutation at residue {data['mutation_resi']} breaks the binding pocket...").
        Keep it under 3 sentences.
        """
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"Gemini could not generate interpretation: {e}"

# --- SIDEBAR (INPUTS) ---
with st.sidebar:
    st.title("LucidDNA 🧬")
    st.caption("Genomic Sentry System")
    
    st.header("1. Upload Data")
    uploaded_file = st.file_uploader("Upload WES/VCF File", type=['txt', 'vcf', 'fastq'])
    
    # Simulation Button
    if st.button("Process Genome", type="primary"):
        with st.status("Initializing Sentry Pipeline...", expanded=True) as status:
            st.write("Aligning to Reference Genome (GRCh38)...")
            time.sleep(1)
            st.write("Filtering 120,000 variants against ClinVar...")
            time.sleep(1)
            st.write("Running ESM-2 Protein Language Model...")
            time.sleep(1)
            status.update(label="Analysis Complete!", state="complete", expanded=False)
        st.session_state['analysis_done'] = True

# --- MAIN DASHBOARD ---
st.title("Preventative Disease Susceptibility Report")

if 'analysis_done' not in st.session_state:
    st.info("Upload a genome file or click 'Process Genome' to start the demo.")
else:
    # 1. METRICS ROW
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total Variants", "124,902")
    c2.metric("Pathogenic Hits", "3", "CRITICAL", delta_color="inverse")
    c3.metric("VUS (Uncertain)", "12", "-2")
    c4.metric("Polygenic Risk", "High", "Type 2 Diabetes")
    
    st.divider()

    # 2. SELECTOR
    st.subheader("Priority Action Items (Triage)")
    selected_gene = st.selectbox("Select a Priority Variant to Inspect:", list(DEMO_DATA.keys()))
    data = DEMO_DATA[selected_gene]

    # 3. GLASS BOX DISPLAY
    col_left, col_right = st.columns([1, 1])

    # Left: Text & AI Analysis
    with col_left:
        st.markdown(f"### 🧬 {selected_gene} Analysis")
        st.markdown(f"**Variant:** `{data['variant']}`")
        st.markdown(f"**Tissue Context:** {data['tissue']}")
        
        st.error(f"**Verdict:** {data['risk']}")
        st.write(f"_{data['description']}_")
        
        st.markdown("#### AI Reliability Score")
        st.progress(abs(data['score'])/20, text=f"ESM-2 Damage Score: {data['score']} (Very High)")

    # Right: 3D Visualization
    with col_right:
        st.markdown(f"### Structural Impact ({data['pdb_id']})")
        
        view, error = render_protein(data['pdb_id'], data['mutation_resi'])
        
        if error:
            st.error(f"Visualization Failed: {error}")
        else:
            # Render the Mol object
            showmol(view, height=400)
            st.caption(f"Real-time render of {selected_gene} structure. Red sphere denotes the specific mutation site.")

    # 4. GEMINI SYNTHESIS
    st.divider()
    with st.expander("Gemini 3.0 Clinical Interpretation (using 2.5 flash rn)", expanded=True):
        # Only fetch if we haven't already for this specific gene (saves API calls)
        if 'gemini_response' not in st.session_state or st.session_state.get('last_gene') != selected_gene:
             with st.spinner("Gemini is analyzing the 3D structure..."):
                 st.session_state['gemini_response'] = get_gemini_interpretation(selected_gene, data)
                 st.session_state['last_gene'] = selected_gene
        
        st.write(st.session_state['gemini_response'])