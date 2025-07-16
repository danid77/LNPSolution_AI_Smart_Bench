from django.urls import path

from . import views

app_name = 'nim'

urlpatterns = [
    path('nimOpenfoldInputPageSample/', views.nimOpenfoldInputPageSample, name='nimOpenfoldInputPageSample'),
    path('nimOpenfoldProcessSample/', views.nimOpenfoldProcessSample, name='nimOpenfoldProcessSample'),
    path('nimOpenfoldOutputPage/', views.nimOpenfoldOutputPage, name='nimOpenfoldOutputPage'),
    path('nimGenmolInputPageSample/', views.nimGenmolInputPageSample, name='nimGenmolInputPageSample'),
    path('nimGenmolProcessSample/', views.nimGenmolProcessSample, name='nimGenmolProcessSample'),
    path('nimGenmolOutputPage/', views.nimGenmolOutputPage, name='nimGenmolOutputPage'),
    path('nimMolminInputPageSample/', views.nimMolminInputPageSample, name='nimMolminInputPageSample'),
    path('nimMolminProcessSample/', views.nimMolminProcessSample, name='nimMolminProcessSample'),
    path('nimMolminOutputPage/', views.nimMolminOutputPage, name='nimMolminOutputPage'),
    path('nimAlphafoldInputPage/', views.nimAlphafoldInputPage, name='nimAlphafoldInputPage'),
    path('nimAlphafoldProcess/', views.nimAlphafoldProcess, name='nimAlphafoldProcess'),
    path('nimAlphafoldOutputPage/', views.nimAlphafoldOutputPage, name='nimAlphafoldOutputPage'),
    path('nimAlphafoldMultiInputPage/', views.nimAlphafoldMultiInputPage, name='nimAlphafoldMultiInputPage'),
    path('nimAlphafoldMultiProcess/', views.nimAlphafoldMultiProcess, name='nimAlphafoldMultiProcess'),
    path('nimAlphafoldMultiOutputPage/', views.nimAlphafoldMultiOutputPage, name='nimAlphafoldMultiOutputPage'),
    path('nimEsmfoldInputPage/', views.nimEsmfoldInputPage, name='nimEsmfoldInputPage'),
    path('nimEsmfoldProcess/', views.nimEsmfoldProcess, name='nimEsmfoldProcess'),
    path('nimEsmOutputPage/', views.nimEsmOutputPage, name='nimEsmOutputPage'),
    path('nimRfdiffusionInputPage/', views.nimRfdiffusionInputPage, name='nimRfdiffusionInputPage'),
    path('nimRfdiffusionProcess/', views.nimRfdiffusionProcess, name='nimRfdiffusionProcess'),
    path('nimRfdiffusionOutputPage/', views.nimRfdiffusionOutputPage, name='nimRfdiffusionOutputPage'),
    path('nimDiffdockInputPage/', views.nimDiffdockInputPage, name='nimDiffdockInputPage'),
    path('nimDiffdockProcess/', views.nimDiffdockProcess, name='nimDiffdockProcess'),
    path('nimDiffdockOutputPage/', views.nimDiffdockOutputPage, name='nimDiffdockOutputPage'),
    path('nimProteinmpnnInputPage/', views.nimProteinmpnnInputPage, name='nimProteinmpnnInputPage'),
    path('nimProteinmpnnProcess/', views.nimProteinmpnnProcess, name='nimProteinmpnnProcess'),
    path('nimProteinmpnnOutputPage/', views.nimProteinmpnnOutputPage, name='nimProteinmpnnOutputPage'),
]