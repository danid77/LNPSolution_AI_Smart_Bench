from django.urls import path

from . import views

app_name = 'plip'

urlpatterns = [
    path('plipInputPage/', views.plipInputPage, name='plipInputPage'),
    path('plipProcess/', views.plipProcess, name='plipProcess'),
    path('plipResultPage/', views.plipResultPage, name='plipResultPage'),
    path('plipResultImages/', views.plipResultImages, name='plipResultImages'),
    path('plipResultXml/', views.plipResultXml, name='plipResultXml'),
]