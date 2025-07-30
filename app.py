import streamlit as st
import pandas as pd
import vcfpy
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import tempfile

# Load ClinVar-style annotations from sample TSV
@st.cache_data
def load_clinvar_annotations(file_path="variant_summary_sample.txt"):
    df = pd.read_csv(file_path, sep='\t', low_memory=False)
    df = df[df['Assembly'] == 'GRCh38']
    db = df.set_index('RS# (dbSNP)')[['GeneSymbol', 'ClinicalSignificance', 'PhenotypeList']].to_dict(orient='index')
    return db

# Annotate each VCF record using the local annotation DB
def annotate_variant(record, annotation_db):
    annotations = []
    variant_id = record.ID
    if not variant_id or variant_id not in annotation_db:
        return annotations
    ann = annotation_db[variant_id]
    for alt in record.ALT:
        annotations.append({
            "CHROM": record.CHROM,
            "POS": record.POS,
            "ID": variant_id,
            "REF": record.REF,
            "ALT": str(alt.value),
            "Gene": ann["GeneSymbol"],
            "Impact": ann["ClinicalSignificance"],
            "Condition": ann["PhenotypeList"]
        })
    return annotations

# Generate a clinical-style PDF report
def generate_pdf_report(annotations, output_pdf):
    c = canvas.Canvas(output_pdf, pagesize=letter)
    width, height = letter
    y = height - 50
    c.setFont("Helvetica-Bold", 14)
    c.drawString(50, y, "Genomic Variant Report")
    y -= 30
    c.setFont("Helvetica", 10)
    for var in annotations:
        line = f"{var['Gene']} | {var['ID']} | {var['Impact']} | {var['Condition']} | {var['CHROM']}:{var['POS']} {var['REF']}>{var['ALT']}"
        c.drawString(50, y, line)
        y -= 20
        if y < 50:
            c.showPage()
            y = height - 50
    c.save()

# Streamlit Web UI
st.title("🧬 Genomic Report Generator")
st.markdown("Upload a **VCF file** to automatically generate a clinical variant PDF report.")

uploaded_file = st.file_uploader("Upload VCF file", type=["vcf"])

if uploaded_file:
    st.success("✅ File uploaded. Parsing...")
    annotation_db = load_clinvar_annotations()

    with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as tmp_vcf:
        tmp_vcf.write(uploaded_file.read())
        tmp_vcf.flush()
        reader = vcfpy.Reader.from_path(tmp_vcf.name)
        annotations = []
        for record in reader:
            ann = annotate_variant(record, annotation_db)
            annotations.extend(ann)

    if annotations:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp_pdf:
            generate_pdf_report(annotations, tmp_pdf.name)
            st.success("📄 PDF report generated.")
            with open(tmp_pdf.name, "rb") as f:
                st.download_button("📥 Download PDF Report", f.read(), file_name="genomic_report.pdf")
    else:
        st.warning("⚠️ No known variants matched the database.")
