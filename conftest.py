"""Put the src/ layout on sys.path so tests can import mechanistic_inference without an install."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
