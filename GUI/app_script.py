from flask import Flask, render_template, request, send_file, session, redirect, url_for
import os
import io
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))  # Add the src directory to the path
from main import generate_assembly_library  # Import the function from main.py

app = Flask(__name__)
app.secret_key = os.urandom(24)

@app.route('/')
def home():
    return render_template('index.html')  # Render the home page with the form

@app.route('/optimize', methods=['POST'])
def optimize():
    # Get form inputs and files
    gene_file = request.files.get('gene_file')
    mutations_file = request.files.get('mutations_file')
    enzyme_name = request.form.get('enzyme_name', 'BbsI')  # Default enzyme
   # min_oligo_size = int(request.form.get('min_oligo_size', 20))
    min_oligo_size = request.form.get('min_oligo_size')
    min_oligo_size = int(min_oligo_size) if min_oligo_size else 20
    max_oligo_size = int(request.form.get('max_oligo_size', 100))

    if not gene_file or not mutations_file:
        return "Error: Both gene file and mutations file are required", 400

    try:
        # Save files temporarily
        temp_dir = './temp'
        os.makedirs(temp_dir, exist_ok=True)
        gene_file_path = os.path.join(temp_dir, gene_file.filename)
        mutations_file_path = os.path.join(temp_dir, mutations_file.filename)
        gene_file.save(gene_file_path)
        mutations_file.save(mutations_file_path)

        # Call the function from main.py
        final_df = generate_assembly_library(
            gene_file_path, mutations_file_path, enzyme_name, min_oligo_size, max_oligo_size
        )

        # Convert the result to CSV for download
        output = io.StringIO()
        final_df.to_csv(output, index=False)
        output.seek(0)
        session['result_csv'] = output.getvalue()  # Save result to session
        return redirect(url_for('result'))
    except Exception as e:
        return f"An error occurred: {str(e)}", 500

@app.route('/gg_assembly')
def gg_assembly():
    return render_template('gg_assembly.html')

@app.route('/result')
def result():
    result_csv = session.get('result_csv', None)
    if not result_csv:
        return redirect(url_for('home'))
    return render_template('result.html')  # Display the result page

@app.route('/download')
def download():
    result_csv = session.get('result_csv', None)
    if not result_csv:
        return redirect(url_for('home'))

    # Send the CSV file as a downloadable response
    output = io.StringIO(result_csv)
    output.seek(0)
    return send_file(output, mimetype='text/csv', as_attachment=True, attachment_filename='result.csv')

if __name__ == '__main__':
    app.run(debug=True)


