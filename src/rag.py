# src/rag.py
import os
from pinecone import Pinecone
from sentence_transformers import SentenceTransformer

def get_rag_context(query):
    # Load keys from environment variables (NOT hardcoded)
    pc = Pinecone(api_key=os.getenv("PINECONE_API_KEY"))
    index = pc.Index("medical-genomics")
    model = SentenceTransformer('all-MiniLM-L6-v2')
    
    # Create embedding and query
    xq = model.encode(query).tolist()
    xc = index.query(vector=xq, top_k=3, include_metadata=True)
    
    return [match['metadata']['text'] for match in xc['matches']]