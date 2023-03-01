from django.db import models

# Create your models here.
from django.db import models
from Bio.Seq import Seq

class Gene(models.Model):
    name = models.CharField(max_length=100)
    sequence = models.TextField()
    created = models.DateTimeField(auto_now_add=True)


    def __str__(self):
        return self.name

    def create_transcript(self, name, sequence):
        transcript = Transcript.objects.create(name=name, sequence=sequence, gene=self)
        return transcript


class Transcript(models.Model):
    name = models.CharField(max_length=100)
    sequence = models.TextField()
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    created = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name
