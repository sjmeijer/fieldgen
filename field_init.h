#ifndef _FIELD_INIT_H
#define _FIELD_INIT_H

/*returns # contacts */
int geometry_init(char *geometry_fn);
enum point_types{OUTSIDE, INSIDE, CONTACT_0, CONTACT_VB};

int init_ev_calc(float ***v[2], float Vb);
int init_wp_calc(float ***v[2], int cno);

float zlen_detector(void);

int geometry_finalize(void);

#endif /*#ifdef _FIELD_INIT_H*/
