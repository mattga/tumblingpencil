import java.util.Timer;

import javax.swing.JFrame;
import javax.vecmath.Vector3d;

public class TumblingPencilSimulation extends JFrame {

	private static final long serialVersionUID = 1L;
	public static Timer timer;
	public static double t, t1, t2, m, m1, m2, r, h1, h2, h, l, Vcylinder, Vcone;
	public static Vector3d p, p1, v, v1i, v1f, omega, omega1i, omega1f, Omega, g;
	public static double [][] Q, Qinv; 
	public static Quaternion q, q1;

	public TumblingPencilSimulation() {
		setUndecorated(true);
		setSize(300, 300);
//		setVisible(true);
		setLocationRelativeTo(null);
	}

	public static void main (String [] args) {

		/**** INPUTS ****/
		// Pencil's Geometry
		t1 = 0.4; t2 = 0.8;
		m = 1.0;
		r = 0.5;
		h1 = 3.0; h2 = 0.5;
		h = (6*h1*h1 + 12*h1*h2 + 3*h2*h2) / (12*h1 + 4*h2);
		l = h1/2 + h2 - h;
		Vcylinder = Math.PI*r*r*h1;
		Vcone = Math.PI*r*r*h2/3;
		m1 = (Vcylinder/(Vcylinder+Vcone))*m; m2 = (Vcone/(Vcylinder+Vcone))*m;

		// Pencil's Velocities
		p1 = new Vector3d(h2*Math.cos(Math.PI/3), 0, h2*Math.sin(Math.PI/3));
		v1i = new Vector3d(-5*Math.cos(Math.PI/6), 0, -5*Math.sin(Math.PI/6));
		omega1i = new Vector3d(1.0, 5, 0.5);
		v1f = new Vector3d(-1.80954, -0.546988, 1.2076);
		omega1f = new Vector3d(0.09957, -0.04174, 0.5);
		g = new Vector3d(0,0,-9.8);

		// Pencil's Moment of Inertia & Rotation Quaternion
		Q = new double[3][3];
		Q[0][0] = Q[1][1] = ((m)/((h1+h2)/3))*(h1*(((3*r*r + h1*h1)/12) + l*l)) + ((h2/3)*(((3/5)*((r*r)/4 + h2*h2)) + h*h));
		Q[2][2] = (.5*m1 + .3*m2)*r*r;
		Qinv = new double[3][3];
		Qinv[0][0] = Qinv[1][1] = 1/Q[0][0];
		Qinv[2][2] = 1/Q[2][2];
		q1 = new Quaternion(Math.cos(Math.PI/12), 0, Math.sin(Math.PI/12), 0);

		Tumble(p1,q1,v1i,omega1i,t1-t2);
		Tumble(p1,q1,v1f,omega1f,t2-t1);
	}

	public static void Tumble (Vector3d _p, Quaternion _q, Vector3d _v, Vector3d _omega, double T) {

		double deltaPhi;
		p = new Vector3d();
		v = new Vector3d();
		
		t = t1;
		omega = _omega;
		q = _q;
		final double delta_t = 0.00001;
		int count = 0;

		while (Math.abs(t-t1) < Math.abs(T)) {
			double t_t1 = t-t1;
			
			if (T < 0) { // Change signs since we're counting down from t=.4 to t=0
				p.x = _p.x - _v.x*t_t1;
				p.y = _p.y - _v.y*t_t1;
				p.z = _p.z - _v.z*t_t1 - .5*g.z*t_t1*t_t1;
			} else {
				p.x = _p.x + _v.x*t_t1;
				p.y = _p.y + _v.y*t_t1;
				p.z = _p.z + _v.z*t_t1 + .5*g.z*t_t1*t_t1;
			}
			
			Omega = q.Lq(omega);
			Vector3d OmegaHat = new Vector3d(Omega);
			OmegaHat.normalize();

			if (T > 0) {
				deltaPhi = Omega.length()*delta_t;
			} else {
				deltaPhi = -1*Omega.length()*delta_t;			
			}

			Quaternion _r = new Quaternion(Math.cos(deltaPhi/2), Math.sin(deltaPhi/2)*OmegaHat.x, Math.sin(deltaPhi/2)*OmegaHat.y, Math.sin(deltaPhi/2)*OmegaHat.z);
			q = Quaternion.multiply(_r, q);

			double [][] omegaArray = {{omega.x},{omega.y},{omega.z}};
			double [][] Qw = MatrixMult(Q, omegaArray);
			Vector3d wxQw = new Vector3d();
			wxQw.cross(omega, new Vector3d(Qw[0][0],Qw[1][0],Qw[2][0]));
			double [][] wxQwArray = {{wxQw.x},{wxQw.y},{wxQw.z}};
			double [][] omegaStep = MatrixMult(Qinv, wxQwArray);
			Vector3d _omegaStep = new Vector3d(omegaStep[0][0],omegaStep[1][0],omegaStep[2][0]);
			_omegaStep.scale(delta_t);
			if (T > 0) {
				omega.sub(_omegaStep);
			} else {
				omega.add(_omegaStep);
			}
			
			if (count++ % 10000 == 0) {
				v.x = _v.x;
				v.y = _v.y;
				if (T < 0) {
					v.z = _v.z - g.z*t_t1;
				} else {
					v.z = _v.z + g.z*t_t1;
				}
				
				System.out.printf("t=%.1f\n", ( T > 0 ? t : t2-t));
				System.out.println("\tp=" + p);
				System.out.println("\tv=" + v);
				System.out.println("\tOmega=" + omega);
			}
			
			t = t + delta_t;
		}
	}

	public static class Quaternion {
		public double q0, q1, q2, q3;

		public Quaternion() {
		}

		public Quaternion (double q0, double q1, double q2, double q3) {
			this.q0 = q0;
			this.q1 = q1;
			this.q2 = q2;
			this.q3 = q3;
		}

		public Vector3d vector() {
			return new Vector3d(q1,q2,q3);
		}

		public Vector3d Lq(Vector3d v) {
			Vector3d result = new Vector3d();

			Quaternion qcc = new Quaternion();
			qcc.q0 = this.q0;
			qcc.q1 = -1 * this.q1;
			qcc.q2 = -1 * this.q2;
			qcc.q3 = -1 * this.q3;
			double s1 = this.q0*this.q0 - Quaternion.multiply(this,qcc).q0;
			double s2 = 2*this.vector().dot(v);
			double s3 = 2*this.q0;
			Vector3d qxv = new Vector3d();
			qxv.cross(this.vector(), v);

			result.x = s1*v.x + s2*this.q1 + s3*qxv.x;
			result.y = s1*v.y + s2*this.q2 + s3*qxv.y;
			result.z = s1*v.z + s2*this.q3 + s3*qxv.z;

			return result;
		}

		public static Quaternion multiply(Quaternion p, Quaternion q) {
			Quaternion r = new Quaternion();
			Vector3d pxq = new Vector3d();
			pxq.cross(p.vector(), q.vector());

			r.q0 = p.q0*q.q0 - p.vector().dot(q.vector());
			r.q1 = p.q0*q.q1 + q.q0*p.q1 + pxq.x;
			r.q2 = p.q0*q.q2 + q.q0*p.q2 + pxq.y;
			r.q3 = p.q0*q.q3 + q.q0*p.q3 + pxq.z;

			return r;
		}

		public String toString() {
			return String.format("%.3f + %.3fi + %.3fj + %.3fk", q0, q1, q2, q3); 
		}
	}

	public static double[][] MatrixMult(double[][] A, double[][] B) {

		int aRows = A.length;
		int aColumns = A[0].length;
		int bRows = B.length;
		int bColumns = B[0].length;

		if (aColumns != bRows) {
			throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
		}

		double[][] C = new double[aRows][bColumns];

		for (int i = 0; i < aRows; i++) { // aRow
			for (int j = 0; j < bColumns; j++) { // bColumn
				for (int k = 0; k < aColumns; k++) { // aColumn
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}

		return C;
	}
}
