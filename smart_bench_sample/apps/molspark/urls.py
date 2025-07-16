from django.urls import path

from . import views

app_name = 'molspark'

urlpatterns = [
    path('', views.molsparkInputPage, name='molsparkInputPage'),
    path('out/', views.molsparkOutputPage, name='molsparkOutputPage'),
    path('molsparkCal/', views.molsparkCal, name='molsparkCal'),
    path('molsparkTable/', views.molsparkTable, name='molsparkTable'),
]