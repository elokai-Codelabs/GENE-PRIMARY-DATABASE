from django.urls import path
from . import views


urlpatterns = [
    path('', views.gene_list, name='gene_list'),
    path('gene/<int:pk>/', views.gene_detail, name='gene_detail'),
    path('search/', views.gene_search, name='gene_search'),
    path('upload/', views.upload_dataset, name='upload'),
    path('add/', views.add_gene, name='add'),
    
]
