from django.urls import path
from . import views
from .views import download_gene_as_pdf




urlpatterns = [
    path('login/', views.loginUser, name='login'),
    path('logout/', views.logoutUser, name='logout'),
    path('', views.index, name='gene_list'),
    path('gene/<int:pk>/', views.gene_detail, name='gene_detail'),
    path('search/', views.gene_search, name='gene_search'),
    path('upload/', views.upload_dataset, name='upload'),
    path('add/', views.add_gene, name='add'),
    path('download_gene/<int:gene_id>/', views.download_gene, name='download_gene'),
    path('gene/<int:gene_id>/download-pdf/', download_gene_as_pdf, name='download_gene_as_pdf'),




    
]
