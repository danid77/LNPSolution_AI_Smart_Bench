from django.urls import path

from . import views

app_name = 'screening'

urlpatterns = [
    path('complexInputPage/', views.complexInputPage, name='complexInputPage'),
    path('complexInputPageSample/', views.complexInputPageSample, name='complexInputPageSample'),
    path('complexOutputPage/', views.complexOutputPage, name='complexOutputPage'),
    # path('screeningInputPage/', views.screeningInputPage, name='screeningInputPage'),
    # path('screeningOutputPage/', views.screeningOutputPage, name='screeningOutputPage'),
    path('generateMolecules/', views.generateMolecules, name='generateMolecules'),
    path('generateMoleculesSample/', views.generateMoleculesSample, name='generateMoleculesSample'),
    path('showMolecules/', views.showMolecules, name='showMolecules'),
]