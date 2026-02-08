import sys
import subprocess
import streamlit as st
import requests
import time
import google.generativeai as genai
from stmol import showmol
import py3Dmol
import pandas as pd
import random
import os

# --- 1. AUTO-INSTALL WINDOWS-FRIENDLY VCF PARSER ---
try:
    import vcfpy
except ImportError:
    st.warning("Installing VCF parser (vcfpy)...")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'vcfpy'])
    import vcfpy

# --- CONFIGURATION ---
st.set_page_config(page_title="LucidDNA Sentry", layout="wide", page_icon="🧬")

# Load API Key
api_key = st.secrets.get("GEMINI_API_KEY", None)
if api_key:
    genai.configure(api_key=api_key)

# --- CONSTANTS: FILE PATHS ---
# These now point specifically to the "data" folder
VCF_PATH = "data/starling_noduprel_qual_miss_filt.recode.vcf"
META_PATH = "data/Metadata_NZ_AU_UK_BE_ReplicatesSibRemoved2.csv"

# --- 2. REAL GENOMIC ANALYSIS ENGINE ---
def analyze_sample_genome(vcf_path, sample_id):
    """
    Parses the VCF to calculate REAL metrics for the specific bird.
    """
    stats = {
        "total_variants": 0,
        "heterozygous": 0,
        "homozygous": 0,
        "max_quality": 0,
        "top_variants": []
    }
    candidate_variants = []

    try:
        if not os.path.exists(vcf_path):
            st.error(f"Analysis Failed: VCF not found at {vcf_path}")
            return None

        reader = vcfpy.Reader.from_path(vcf_path)
        
        # Check if sample exists in VCF header
        if sample_id not in reader.header.samples.names:
            st.error(f"Sample ID '{sample_id}' not found in VCF header.")
            return None

        # Iterate through VCF records
        for record in reader:
            call = record.call_for_sample.get(sample_id)
            
            if call and call.is_variant:
                stats["total_variants"] += 1
                
                if call.is_het:
                    stats["heterozygous"] += 1
                else:
                    stats["homozygous"] += 1
                
                # Capture Quality
                qual = record.QUAL if record.QUAL else 0
                if qual > stats["max_quality"]:
                    stats["max_quality"] = qual
                
                # Save high-quality variants for dropdown
                if qual > 20: 
                    candidate_variants.append({
                        "id": f"{record.CHROM}:{record.POS}",
                        "chrom": record.CHROM,
                        "pos": record.POS,
                        "ref": record.REF,
                        "alt": record.ALT[0].value,
                        "qual": qual,
                        "description": f"Type: {'SNP' if len(record.REF)==1 else 'Indel'}"
                    })
        
        # Sort by Quality and pick top 3
        candidate_variants.sort(key=lambda x: x['qual'], reverse=True)
        stats["top_variants"] = candidate_variants[:3]
        
        # Calculate Ratio
        if stats["total_variants"] > 0:
            stats["het_ratio"] = stats["heterozygous"] / stats["total_variants"]
        else:
            stats["het_ratio"] = 0
            
    except Exception as e:
        st.error(f"Analysis Error: {e}")
        return None

    return stats

# --- 3. DATA LOADING (SIDEBAR) ---
@st.cache_data
def load_data_safe():
    try:
        # Debugging: Check if files exist
        if not os.path.exists(META_PATH):
            st.error(f"❌ Metadata not found at: {META_PATH}")
            return None, None
        if not os.path.exists(VCF_PATH):
            st.error(f"❌ VCF not found at: {VCF_PATH}")
            return None, None

        df = pd.read_csv(META_PATH)
        
        # Fast text parsing for sample list (Avoiding huge VCF load time)
        samples = []
        with open(VCF_PATH, "r") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    samples = line.strip().split("\t")[9:]
                    break
        
        return df, samples
    except Exception as e:
        st.error(f"Data Load Error: {e}")
        return None, None

metadata_df, vcf_samples = load_data_safe()

# --- 4. CALLBACK FOR ROBUST RANDOMIZATION ---
def randomize_callback():
    """Forces the update before the page redraws."""
    if metadata_df is not None and vcf_samples:
        selected_id = random.choice(vcf_samples)
        
        # Try to find matching metadata
        # Strategy: exact match OR match without .sorted.bam extension
        row = metadata_df[metadata_df['id'] == selected_id]
        if row.empty:
            clean_id = selected_id.replace(".sorted.bam", "")
            row = metadata_df[metadata_df['id'].str.contains(clean_id, regex=False)]
        
        if not row.empty:
            data = row.iloc[0]
            st.session_state['selected_id'] = selected_id # Keep exact VCF ID for analysis
            st.session_state['origin'] = f"{data['pop']} ({data['Con']})"
            st.session_state['lat_lon'] = [data['lat'], data['lon']]
            st.session_state['analysis_complete'] = False # Reset dashboard
        else:
            # Fallback if metadata is missing for this specific bird
            st.session_state['selected_id'] = selected_id
            st.session_state['origin'] = "Unknown Location"
            st.session_state['lat_lon'] = None
            st.session_state['analysis_complete'] = False

# --- HELPER FUNCTIONS (Gemini & 3D) ---
def render_protein(pdb_id="4HHB", resi=1):
    try:
        url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            view = py3Dmol.view(width=500, height=400)
            view.addModel(response.text, 'pdb')
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.addStyle({'resi': str(resi)}, {'sphere': {'color': 'red'}})
            view.zoomTo()
            return view
    except:
        pass
    return None

def get_gemini_analysis(variant, metrics):
    if not api_key: return "API Key Missing."
    prompt = f"""
    Analyze this Starling genome variant (Conservation Genetics Context).
    
    Variant: {variant['chrom']} at {variant['pos']} ({variant['ref']} -> {variant['alt']})
    Quality Score: {variant['qual']}
    Sample Heterozygosity: {metrics['het_ratio']:.3f}
    
    Provide a 2-sentence assessment of:
    1. The reliability of this call (based on QUAL).
    2. Implications for genetic diversity (is the bird inbred?).
    """
    try:
        model = genai.GenerativeModel('gemini-2.5-flash')
        return model.generate_content(prompt).text
    except Exception as e:
        return f"AI Error: {e}"

# --- SIDEBAR ---
with st.sidebar:
    st.title("LucidDNA 🧬")
    st.caption("Sentry Mode: Active")
    
    st.header("1. Sample Selection")
    if vcf_samples:
        st.write(f"Database: {len(vcf_samples)} Individuals")
        st.button("🎲 Randomize Subject", on_click=randomize_callback)
    else:
        st.error("Database connection failed.")

    if 'selected_id' in st.session_state:
        st.success(f"**Target:** {st.session_state['selected_id']}")
        st.info(f"**Origin:** {st.session_state['origin']}")
    
    st.header("2. Analysis")
    if st.button("🚀 Run Sentry Pipeline", type="primary", disabled='selected_id' not in st.session_state):
        with st.status("Processing Genomic Data...", expanded=True) as status:
            st.write(f"📂 Mounting {VCF_PATH}...")
            time.sleep(0.5)
            
            # CALLING THE ANALYSIS WITH CORRECT PATH
            results = analyze_sample_genome(VCF_PATH, st.session_state['selected_id'])
            
            if results:
                st.session_state['genome_stats'] = results
                status.update(label="Sequencing Complete", state="complete", expanded=False)
                st.session_state['analysis_complete'] = True
            else:
                status.update(label="Analysis Failed", state="error")

# --- MAIN DASHBOARD ---
st.title("Genomic Integrity Report")

if not st.session_state.get('analysis_complete'):
    st.info("👈 Select a subject in the sidebar and run the pipeline.")
else:
    stats = st.session_state['genome_stats']
    
    # 1. METRICS
    risk_score = stats['het_ratio']
    verdict = "HIGH RISK (Inbred)" if risk_score < 0.2 else "STABLE"
    
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Loci Analyzed", f"{stats['total_variants']:,}")
    c2.metric("Heterozygosity", f"{stats['het_ratio']:.3f}")
    c3.metric("Max Quality", f"{stats['max_quality']:.1f}")
    c4.metric("Verdict", verdict, delta=f"{risk_score:.2f}")
    
    # 3. DETAILED VARIANT ANALYSIS
    st.divider()
    st.subheader("⚠️ High-Impact Variants Identified")
    
    if stats['top_variants']:
        # Create friendly labels for dropdown
        opts = {f"{v['chrom']}:{v['pos']} ({v['ref']}->{v['alt']})" : v for v in stats['top_variants']}
        sel_label = st.selectbox("Select Variant for AI Interpretation:", list(opts.keys()))
        sel_data = opts[sel_label]
        
        c_left, c_right = st.columns([1,1])
        with c_left:
            st.markdown(f"**Locus:** {sel_data['chrom']} - {sel_data['pos']}")
            st.markdown(f"**Confidence:** {sel_data['qual']:.1f}")
            st.markdown("#### Gemini Analysis")
            with st.spinner("Analyzing conservation impact..."):
                st.write(get_gemini_analysis(sel_data, stats))
        
        with c_right:
            view = render_protein("4HHB", resi=15) # Placeholder structure
            if view:
                showmol(view, height=350)
                st.caption("Structural Homolog Visualization")
    else:
        st.warning("No high-quality variants found in this sample.")