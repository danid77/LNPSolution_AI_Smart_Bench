from django.urls import path

from . import views

app_name = 'visualzing'

urlpatterns = [
    path('', views.visualzingInputPage, name='visualzingInputPage'),
    path('visualzingProcess/', views.visualzingProcess, name='visualzingProcess'),
    path('visualzingResultPage/', views.visualzingResultPage, name='visualzingResultPage'),
    path('visualzingPdb/', views.visualzingPdb, name='visualzingPdb'),
    path('visualzingResult2DImage/', views.visualzingResult2DImage, name='visualzingResult2DImage'),
    path('visualzingPlipResultImages/', views.visualzingPlipResultImages, name='visualzingPlipResultImages'),
    # path('visualzingPlipResultPdbs/', views.visualzingPlipResultPdbs, name='visualzingPlipResultPdbs'),
    path('visualzingPlipLigand/', views.visualzingPlipLigand, name='visualzingPlipLigand'),
    path('visualzingPlipScore/', views.visualzingPlipScore, name='visualzingPlipScore'),
    path('visualzingGninaScore/', views.visualzingGninaScore, name='visualzingGninaScore'),
]