import re
from .models import Gene
from Bio.Seq import Seq

from django.http import HttpResponse
# this is the main views.py
from django.shortcuts import render,redirect, get_object_or_404
from django.http import HttpResponse
from django.contrib import messages
# fasta 
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO
from .models import Gene, Transcript

from django.contrib.auth import login, authenticate, logout
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required

from django.shortcuts import render, get_object_or_404
from django.utils.translation import gettext_lazy as _
# regex 

from django.template.loader import get_template
from django.shortcuts import get_object_or_404
from xhtml2pdf import pisa
from io import BytesIO


# ============USER AUTHENTICATION , LOGIN AND LOGOUT==========
# ============================================================
def loginUser(request):
   
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']


        user = authenticate(request, username=username, password=password)

        if user is not None:
            login(request, user)
            messages.success(request,'User succesfully logged in')
            return redirect('gene_list')
        else:
            messages.error(request,'Invalid username or password')
    context ={}
    return render(request,'app/login.html', context)   


@login_required(login_url='login')
def logoutUser(request):
    logout(request)
    return redirect('login')


# @login_required(login_url='login')
def index(request):
    genes = Gene.objects.all()
    genes_count = Gene.objects.count()
    return render(request, 'app/index.html', {'genes': genes,'genes_count':genes_count})


@login_required(login_url='login')
def add_gene(request):
    if request.method == 'POST':
        name = request.POST['name']
        sequence = request.POST['sequence']

        # Validate sequence input
        pattern = re.compile('^[ATCG]*$')
        if not pattern.match(sequence):
            messages.error(request, 'Invalid sequence. Only A, T, C, and G are allowed.')
            return render(request, 'app/add_gene.html')

        Gene.objects.create(name=name, sequence=sequence)
        messages.success(request, 'Gene added successfully.')
        return redirect('gene_list')
    else:
        return render(request, 'app/add_gene.html')
    




def gene_detail(request, pk):
    gene = get_object_or_404(Gene, pk=pk)
    trans_gene = Seq(str(gene))
    gene_sequence =  str(gene.sequence)
    mRNA_sequence =Seq(gene_sequence).complement_rna()
    protein_sequence = Seq(gene_sequence).complement_rna().translate()

    context = {'gene': gene, 'transcripts':trans_gene,'mRNA_sequence':mRNA_sequence,'protein_sequence':protein_sequence}
    return render(request, 'app/gene_detail.html', context)

    
def gene_search(request):
    gene_name = request.GET.get('gene_name')
    genes = Gene.objects.filter(name__icontains=gene_name)
    return render(request, 'app/gene_list.html', {'genes': genes})

@login_required(login_url='login')
def upload_dataset(request):
    if request.method == 'POST' and request.FILES['fasta_file']:
        # Get the uploaded file from the form
        uploaded_file = request.FILES['fasta_file']

        # Save the file to disk using Django's file storage API
        fs = FileSystemStorage()
        filename = fs.save(uploaded_file.name, uploaded_file)

        # Open the file using BioPython's SeqIO.parse() function
        filepath = fs.path(filename)
        with open(filepath) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Create a Gene object for each record in the file
                gene = Gene.objects.create(name=record.id, sequence=str(record.seq))

                # Create a Transcript object for each record in the file
                transcript = Transcript.objects.create(name=record.id, sequence=str(record.seq), gene=gene)

        # Delete the uploaded file from disk
        fs.delete(filename)

        # Redirect the user back to the home page
        return redirect('gene_list')

    # If the form was not submitted or the file was not provided, render the upload form
    return render(request, 'app/upload.html')


def download_gene(request, gene_id):
    # Retrieve the gene from the database
    gene = Gene.objects.get(id=gene_id)

    # Format the gene as FASTA
    fasta_data = ">{}\n{}\n".format(gene.name, gene.sequence)

    # Return the gene as a file download
    response = HttpResponse(fasta_data, content_type="text/plain")
    response["Content-Disposition"] = "attachment; filename={}.fasta".format(gene.name)
    return response


def download_gene_as_pdf(request, gene_id):
    # Retrieve the gene from the database
    gene = get_object_or_404(Gene, id=gene_id)

    # Format the gene data as a dictionary
    gene_data = {'gene_name': gene.name, 'gene_sequence': gene.sequence}

    # Create a PDF template using the gene data
    template = f"""
    <html>
        <head>
            <style>
                h1 {{
                    color: #333;
                    font-size: 24px;
                    text-align: center;
                    margin-bottom: 24px;
                }}
                p {{
                    font-size: 16px;
                    line-height: 1.5;
                    margin: 0;
                    padding: 0;
                    height: 500px; /* Set a fixed height for the container */
                    overflow: auto; /* Add overflow property to scroll if content exceeds container height */
                    white-space: pre-wrap;
                    word-wrap: break-word;
                }}
            </style>
        </head>
        <body>
            <h1>{gene.name}</h1>
            <p>{gene.sequence}</p>
        </body>
    </html>
    """

    pdf_file = BytesIO()
    pisa.CreatePDF(BytesIO(template.encode("UTF-8")), dest=pdf_file)

    # Return the PDF as a file download
    pdf_file.seek(0)
    response = HttpResponse(pdf_file, content_type='application/pdf')
    response['Content-Disposition'] = 'attachment; filename="{}.pdf"'.format(gene.name)
    return response






