from django.urls import path

from . import views

app_name = 'lnp'

urlpatterns = [
    path('', views.home),
    path('workFlow/', views.workFlow, name='workFlow'),
    path('singleModule/', views.singleModule, name='singleModule'),
    path('proteinHome/', views.proteinHome, name='proteinHome'),
    path('moleculeHome/', views.moleculeHome, name='moleculeHome'),
    path('complexHome/', views.complexHome, name='complexHome'),
    path('peptideHome/', views.peptideHome, name='peptideHome'),
    path('interactionHome/', views.interactionHome, name='interactionHome'),
    path('dimensionalityHome/', views.dimensionalityHome, name='dimensionalityHome'),
    path('test/', views.test, name='test'),
    # path('lnpWork/', views.lnpWork, name='lnpWork'),
    # path('nvidiaWork/', views.nvidiaWork, name='nvidiaWork'),
    # path('mix/', views.mix, name='mix'),
]