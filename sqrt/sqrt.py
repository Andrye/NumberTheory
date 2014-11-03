def mod_exp(a,d,p):
	q = 1
	while d > 0:
		if d % 2 == 1:
			q=(q * a)%p
        	a = (a * a)%p
        	d /= 2
	return q


'''
    Returns a solution 's' to the congruence of a form n² ≡ s (mod p) 
'''

def Tonelli_Shanks(n,p):
	if p == 2:
		return n
	if mod_exp(n, (p-1)/2, p) != 1:
		return "ERROR"
	d = p-1 
	Q, S = d, 0
	while Q%2 == 0:
		S, Q = S+1, Q/2
	z = 2
	while mod_exp(z, d/2, p) != d:
		z = z+1
	c, R, t, M=mod_exp(z, Q, p), mod_exp(n, (Q+1)/2, p), mod_exp(n, Q, p), S
	while t!=1:
		i=1
		tq=(t**2)%p
		while tq!=1:
			tq, i = (tq*tq)%p, i+1
		b = mod_exp(c, 2**(M-i-1), p)
		R = (R*b)%p
		t, c, M = (t*b*b)%p, (b*b)%p, i
	return R

