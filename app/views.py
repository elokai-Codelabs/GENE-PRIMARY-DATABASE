from .models import Gene
from Bio.Seq import Seq
# this is the main views.py
from django.shortcuts import render,redirect, get_object_or_404
from django.http import HttpResponse
from django.contrib import messages
# fasta 
from django.core.files.storage import FileSystemStorage
from Bio import SeqIO
from .models import Gene, Transcript

def gene_list(request):
    genes = Gene.objects.all()
    genes_count = Gene.objects.count()
    return render(request, 'app/index.html', {'genes': genes,'genes_count':genes_count})


def add_gene(request):
    if request.method == 'POST':
        name = request.POST['name']
        sequence = request.POST['sequence']
        Gene.objects.create(name=name, sequence=sequence)
        messages.success(request, 'Gene added successfully.')
        return redirect('gene_list')
    else:
        return render(request, 'app/add_gene.html')

def gene_detail(request, pk):
    gene = get_object_or_404(Gene, pk=pk)
    trans_gene = Seq(str(gene)).complement_rna()
    return render(request, 'app/gene_detail.html', {'gene': gene,'transcripts':trans_gene})

    

def gene_search(request):
    gene_name = request.GET.get('gene_name')
    genes = Gene.objects.filter(name__icontains=gene_name)
    return render(request, 'app/gene_list.html', {'genes': genes})


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



