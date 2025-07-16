"""
URL configuration for config project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include

from apps.lnp import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('lnp/', include('apps.lnp.urls')),
    path('molspark/', include('apps.molspark.urls')),
    path('screening/', include('apps.screening.urls')),
    path('plip/', include('apps.plip.urls')),
    path('visualzing/', include('apps.visualzing.urls')),
    path('umap/', include('apps.umap.urls')),
    path('nim/', include('apps.nim.urls')),
    path('protein/', include('apps.protein.urls')),
]
