#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TAM 101
#define STEPS 3000001
#define real double

int main()
{
    FILE *saida = fopen("hipertermia_heterogeneo_np_0.txt", "w");
    if (saida == NULL)
    {
        printf("Ocorreu um erro na abertura do arquivo de saída.\n");
        return 1;
    }
    int i, j, k, w, s, swap = 0;
    real pb = 1000, Ta = 37, cb = 4200, **p, **c, **Q, A = 1.3*pow(10, 6), r0 = 3.1*pow(10, -3);
    real ***u, uim, uip, ujp, ujm, hx = 0.001, ht = 0.001, **wb, **Calor_nano, **kappa, kip, kim, kjp, kjm;
    char str[40];
    Q = malloc(TAM*sizeof(real*));
    Calor_nano = malloc(TAM*sizeof(real*));
    p = malloc(TAM*sizeof(real*));
    u = malloc(2*sizeof(real**));
    c = malloc(TAM*sizeof(real*));
    wb = malloc(TAM*sizeof(real*));
    kappa = malloc(TAM*sizeof(real*));
    for (i = 0; i < TAM; i++)
    {
        Q[i] = malloc(TAM*sizeof(real));
        Calor_nano[i] = malloc(TAM*sizeof(real));
        p[i] = malloc(TAM*sizeof(real));
        c[i] = malloc(TAM*sizeof(real));
        wb[i] = malloc(TAM*sizeof(real));
        kappa[i] = malloc(TAM*sizeof(real));
    }
    for (i = 0; i < 2; i++)
        u[i] = malloc(TAM*sizeof(real*));
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            u[i][j] = malloc(TAM*sizeof(real));
        }
    }
    //definindo as constantes
    for (i = 0; i < TAM; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            u[swap][i][j] = Ta;
            Calor_nano[i][j] = A*exp(-(pow((i*hx - 0.05), 2) + pow((j*hx - 0.05), 2))/pow(r0, 2));
            if ((i*hx >= 0.04 && i*hx <= 0.06) && (j*hx >= 0.04 && j*hx <= 0.06))
            {
                p[i][j] = 1000;
                c[i][j] = 4200;
                Q[i][j] = 4200;
                wb[i][j] = 1.25*(pow(10, -3));
                kappa[i][j] = 0.55;
            }
              else
            {
                p[i][j] = 1200;
                c[i][j] = 3600;
                Q[i][j] = 420;
                wb[i][j] = 5*(pow(10, -4));
                kappa[i][j] = 0.5;
            }
        }
    }
    for (i = 0; i < TAM; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            fprintf(saida, "%lf ", u[swap][i][j]);
        }
        fprintf(saida, "\n");
    }
    fclose(saida);
    for (k = 0; k < STEPS; k++)
    {
      for (i = 0; i < TAM; i++)
      {
          for (j = 0; j < TAM; j++)
          {
              switch(i)
              {
                  case 0:
                    {
                        uim = u[swap][1][j];
                        uip = u[swap][i+1][j];
                        kim = kappa[i][j];
                        kip = 2*kappa[i+1][j]*kappa[i][j]/(kappa[i][j] + kappa[i+1][j]);
                        break;
                    }
                  case (TAM-1):
                    {
                        uip = u[swap][TAM-2][j];
                        uim = u[swap][i-1][j];
                        kip = kappa[i][j];
                        kim = 2*kappa[i-1][j]*kappa[i][j]/(kappa[i-1][j] + kappa[i][j]);
                        break;
                    }
                  default:
                    {
                        uim = u[swap][i-1][j];
                        uip = u[swap][i+1][j];
                        kip = 2*kappa[i+1][j]*kappa[i][j]/(kappa[i][j] + kappa[i+1][j]);
                        kim = 2*kappa[i-1][j]*kappa[i][j]/(kappa[i-1][j] + kappa[i][j]);
                        break;
                    }
              }
              switch(j)
              {
                  case 0:
                    {
                        ujm = u[swap][i][1];
                        ujp = u[swap][i][j+1];
                        kjm = kappa[i][j];
                        kjp = 2*kappa[i][j+1]*kappa[i][j]/(kappa[i][j] + kappa[i][j+1]);
                        break;
                    }
                  case (TAM-1):
                    {
                        ujp = u[swap][i][TAM-2];
                        ujm = u[swap][i][j-1];
                        kjp = kappa[i][j];
                        kjm = 2*kappa[i][j-1]*kappa[i][j]/(kappa[i][j-1] + kappa[i][j]);
                        break;
                    }
                  default:
                    {
                        ujm = u[swap][i][j-1];
                        ujp = u[swap][i][j+1];
                        kjp = 2*kappa[i][j+1]*kappa[i][j]/(kappa[i][j] + kappa[i][j+1]);
                        kjm = 2*kappa[i][j-1]*kappa[i][j]/(kappa[i][j-1] + kappa[i][j]);
                        break;
                    }
              }
              u[!swap][i][j] = (ht/(p[i][j]*c[i][j]))*((kip*uip + kjp*ujp - (kip+kim+kjp+kjm)*u[swap][i][j] + kim*uim + kjm*ujm)/(pow(hx, 2)) + wb[i][j]*pb*cb*(Ta - u[swap][i][j]) + Q[i][j] + Calor_nano[i][j]) + u[swap][i][j];
          }
      }
      if (k % 300000 == 0)
      {
          sprintf(str, "hipertermia_heterogeneo_np_%d.txt", (k/300000)+1);
          FILE *saida = fopen(str, "w");
      }
      for (w = 0; w < TAM; w++)
      {
          for (s = 0; s < TAM; s++)
          {
              u[swap][w][s] = u[!swap][w][s];
              if (k % 300000 == 0)
              {
                  fprintf(saida, "%lf ", u[swap][w][s]);
              }
          }
          if (k % 300000 == 0)
            fprintf(saida, "\n");
      }
      if (k % 300000 == 0)
      {
          fclose(saida);
      }
      swap = !swap;
    }
    for (i = 0; i < TAM; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            printf("%lf ", u[swap][i][j]);
        }
        printf("\n");
        free(Q[i]);
        free(p[i]);
        free(c[i]);
        free(wb[i]);
        free(Calor_nano[i]);
        free(kappa[i]);
    }
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < TAM; j++)
            free(u[i][j]);
        free(u[i]);
    }
    free(Q);
    free(p);
    free(u);
    free(c);
    free(wb);
    free(Calor_nano);
    free(kappa);
    return 0;
}

