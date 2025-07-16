from django.urls import path

from . import views

app_name = 'protein'

urlpatterns = [
    path('proteinInputPage/', views.proteinInputPage, name='proteinInputPage'),
    path('proteinInputPageSample/', views.proteinInputPageSample, name='proteinInputPageSample'),
    path('proteinProcess/', views.proteinProcess, name='proteinProcess'),
    path('proteinProcessSample/', views.proteinProcessSample, name='proteinProcessSample'),
    path('proteinOutputPage/', views.proteinOutputPage, name='proteinOutputPage'),
    path('boltzInputPage/', views.boltzInputPage, name='boltzInputPage'),
    path('boltzProcess/', views.boltzProcess, name='boltzProcess'),
    path('boltzOutputPage/', views.boltzOutputPage, name='boltzOutputPage'),
    path('chaiInputPage/', views.chaiInputPage, name='chaiInputPage'),
    path('chaiProcess/', views.chaiProcess, name='chaiProcess'),
    path('chaiOutputPage/', views.chaiOutputPage, name='chaiOutputPage'),
    path('boltzComplexInputPage/', views.boltzComplexInputPage, name='boltzComplexInputPage'),
    path('boltzComplexProcess/', views.boltzComplexProcess, name='boltzComplexProcess'),
    path('boltzComplexOutputPage/', views.boltzComplexOutputPage, name='boltzComplexOutputPage'),
    path('chaiComplexInputPage/', views.chaiComplexInputPage, name='chaiComplexInputPage'),
    path('chaiComplexProcess/', views.chaiComplexProcess, name='chaiComplexProcess'),
    path('chaiComplexOutputPage/', views.chaiComplexOutputPage, name='chaiComplexOutputPage'),
]