// -*- C++ -*-

//���̃v���O�����́A�������ʂ�蓾��ꂽRange�Atheta�Aphi���A���݂��m�F����Ă��邷�ׂĂ̗��q�ɂ�����track�𓖂Ă͂߂āA�e���q�̎��ʂ��v�Z����B
//���ꂪ�덷�̂R�{�ȓ��Ő��藧���̂������o�����̂ł���B�\������ۂ́A�ǂ̗��q�𓖂Ă͂߂āA�v�Z���������AZ,A,S��p���āA�\�����Ă��܂��B
//�o�͂���result.txt�ɂ́A��������\�L���āA�ǂ̂悤�Ȕ��������肳�ꂽ�������₷�����Ă���܂��B

#include <stdio.h>
#include <math.h>

#include "charged.hh"
#include "n_charged.hh"
#include "P.hh"
#include "KE.hh"

//______________________________________________________________________________
int
main( int argc, char *argv[] )
{
  int    charged_kind, n_kind;
  double Range1, Range2, Range3, errRange1, errRange2, errRange3;
  double theta1, theta2, theta3, errtheta1, errtheta2, errtheta3;
  double phi1, phi2, phi3, errphi1, errphi2, errphi3;

  //�L�[�{�[�h��������l����͂���B
  printf("\n\n");
  printf("�i�P�j���̃v���O�����́A�L�[�{�[�h���ł����܂ꂽtrack��Range�A�ƁA�ӂ���\n");
  printf("�@�@�@���Ă͂܂����ߒ������������v���O�����ł��B\n");
  printf("�i�Q�j���̃v���O�����ł́A�ȉ��̕���p�^�[�����A���l���͂ɂđI�����A\n");
  printf("      �g�p���邱�Ƃ��ł��܂��B\n");
  printf("�i�R�j���炩���ߌ��߂Ă������A���q�P���痱�q�R�̏��ɑł�����ł��������B\n");
  printf("�i�S�j�v�Z���������܂�����A���̃E�C���h�E�Ɂu�v�Z����܂����I�v�ƕ\�L����܂��B");
  printf("�@�@�@�\�L�����܂ŁA���΂炭���҂����������B\n");
  printf("�i�T�j�ł����ސ��l�́A���p�����œ��͂��Ă��������B\n");
  printf("�i�U�j���l�́A�ł����ނ��Ƃ�enter�L�[�������Ă��������B\n\n\n");

  printf("�͂��߂�,����p�^�[�������肵�܂��B\n");
  printf("�ucharged particle �̐��v�́A�u�P�v�`�u�R�v����͂��Ă��������B\n");
  printf("�u�������q�̐��v�́A�u�O�v���u�P�v����͂��Ă��������B\n\n\n");

  printf("charged particle �̐� = ");
  scanf("%d", &charged_kind);

  printf("\n�������q�̐� = ");
  scanf("%d", &n_kind);



  printf("\n\n��������́ARange�A�ƁA�ӂ̓��͂��s���܂��B\n\n\n");

  printf("Range�P= ");
  scanf("%lf", &Range1);
  printf("��Range�P= ");
  scanf("%lf", &errRange1);
  printf("�ƂP = ");
  scanf("%lf", &theta1);
  printf("�ЃƂP = ");
  scanf("%lf", &errtheta1);
  printf("�ӂP= ");
  scanf("%lf", &phi1);
  printf("�ЃӂP= ");
  scanf("%lf", &errphi1);
  printf("\nRange�Q = ");
  scanf("%lf", &Range2);
  printf("��Range�Q = ");
  scanf("%lf", &errRange2);
  printf("�ƂQ = ");
  scanf("%lf", &theta2);
  printf("�ЃƂQ = ");
  scanf("%lf", &errtheta2);
  printf("�ӂQ= ");
  scanf("%lf", &phi2);
  printf("�ЃӂQ= ");
  scanf("%lf", &errphi2);
  printf("\nRange�R= ");
  scanf("%lf", &Range3);
  printf("��Range�R= ");
  scanf("%lf", &errRange3);
  printf("�ƂR = ");
  scanf("%lf", &theta3);
  printf("�ЃƂR = ");
  scanf("%lf", &errtheta3);
  printf("�ӂR= ");
  scanf("%lf", &phi3);
  printf("�ЃӂR= ");
  scanf("%lf", &errphi3);




  //Range1 = 68.5;
  //Range2 = 37.4;
  //Range3 = 0;
  //errRange1 = 1.6;
  //errRange2 = 1.3;
  //errRange3 = 0;
  //theta1 = 117.9;
  //theta2 = 61.2;
  //theta3 = 0;
  //errtheta1 = 3;
  //errtheta2 = 7;
  //errtheta3 = 0;
  //phi1 = 179.9;
  //phi2 = 0.1;
  //phi3 = 0;
  //errphi1 = 0.2;
  //errphi2 = 7;
  //errphi3 = 0;


  //�����Ă���̂��m�F�B
  printf("\n\n�^����ꂽ�l�����ƂɁA�ȉ��Ɏ��������e�Ōv�Z���܂��B\n\n");


  if((charged_kind==3)&&(n_kind==0))//charged particles ���R�̏ꍇ
    {
      printf("\n\ncharged particles ���R�{�̌v�Z���s���܂��B\n\n");

      double answer1 = three_charged(Range1, Range2, Range3, errRange1, errRange2, errRange3, theta1, theta2, theta3, errtheta1, errtheta2, errtheta3, phi1, phi2, phi3, errphi1, errphi2, errphi3);

      printf("\n\n�v�Z���I�����܂���!");
      getchar();
    }

  else if((charged_kind==2)&&(n_kind==0))//charged partibles ���Q�̏ꍇ
    {
      printf("\n\ncharged particles ���Q�{�̌v�Z���s���܂��B\n\n");

      double answer2 = two_charged(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n�v�Z���I�����܂���!");
      getchar();
    }


  else if((charged_kind==1)&&(n_kind==1))//charged particles ���P�ƒ������q���P�̏ꍇ
    {
      printf("\n\ncharged particles ��1�{�ƁA�������q���P�i�����q�̏ꍇ�A�Q���܂ށj�̌v�Z���s���܂��B\n\n");

      double answer3 = n_one(Range1, errRange1, theta1, errtheta1, phi1, errphi1);

      printf("\n\n�v�Z���I�����܂����I");
      getchar();

    }


  else if((charged_kind==2)&&(n_kind==1))//charged particles ���Q�ƒ������q���P�̏ꍇ
    {
      printf("\n\ncharged particles ���Q�{�ƁA�������q���P�i�����q�̏ꍇ�A�Q���܂ށj�̌v�Z���s���܂��B\n\n");

      double answer5 = n_two(Range1, Range2, errRange1, errRange2, theta1, theta2, errtheta1, errtheta2, phi1, phi2, errphi1, errphi2);

      printf("\n\n�v�Z���I�����܂����I");
      getchar();

    }


  else if((charged_kind==3)&&(n_kind==1))//charged particles ���R�ƒ������q���P�̏ꍇ
    {
      printf("\n\ncharged particles ���R�{�ƁA�������q���P�i�����q�̏ꍇ�A�Q���܂ށj�̌v�Z���s���܂��B\n\n");

      double answer5 = n_three(Range1, Range2, Range3, errRange1, errRange2, errRange3, theta1, theta2, theta3, errtheta1, errtheta2, errtheta3, phi1, phi2, phi3, errphi1, errphi2, errphi3);

      printf("\n\n�v�Z���I�����܂����I");
      getchar();

    }

}
