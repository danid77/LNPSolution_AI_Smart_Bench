from django.urls import path

from . import views

app_name = 'umap'

urlpatterns = [
    path('', views.umapInputPage, name='umapInputPage'),
    path('umapInputPageSample/', views.umapInputPageSample, name='umapInputPageSample'),
    path('generateMolecules/', views.generateMolecules, name='generateMolecules'),
    path('generateMoleculesSample/', views.generateMoleculesSample, name='generateMoleculesSample'),
    path('umapResultPage/', views.umapResultPage, name='umapResultPage'),
    path('generate3dSdf/', views.generate3dSdf, name='generate3dSdf'),
    path('umapResultToolNum/', views.umapResultToolNum, name='umapResultToolNum'),
    path('umapResultImage/', views.umapResultImage, name='umapResultImage'),
    path('umapResultCsv/', views.umapResultCsv, name='umapResultCsv'),
]