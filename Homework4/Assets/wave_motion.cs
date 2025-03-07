using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;


	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{		
		//Step 1:
		//TODO: Compute new_h based on the shallow wave model.
		for(int i = 0; i < size; i ++)
        {
			for(int j = 0; j < size; j ++)
            {
				new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping;
                new_h[i, j] += compute_pressure_with_neumann_boundary(h, i, j, i -1, j);
				new_h[i, j] += compute_pressure_with_neumann_boundary(h, i, j, i + 1, j);
				new_h[i, j] += compute_pressure_with_neumann_boundary(h, i, j, i, j - 1);
				new_h[i, j] += compute_pressure_with_neumann_boundary(h, i, j, i, j + 1);
            }
        }

		//Step 2: Block->Water coupling
		//TODO: for block 1, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).

		//TODO: for block 2, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
	
		//TODO: Diminish vh.

		//TODO: Update new_h by vh.

		//Step 3
		//TODO: old_h <- h; h <- new_h;
		for(int i = 0; i < size; i ++)
        {
			for(int j = 0; j < size; j ++)
            {
				old_h[i, j] = h[i, j];
				h[i, j] = new_h[i, j];
            }
        }

		//Step 4: Water->Block coupling.
		//More TODO here.
	}
	float compute_pressure_with_neumann_boundary(float[,] h, int i , int j, int i_s, int j_s)
    {
		if (i_s >= 0 && i_s < size && j_s >= 0 && j_s < size)
			return (h[i_s, j_s] - h[i, j]) * rate;
		return 0;
    }
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

		//TODO: Load X.y into h.
		for(int v = 0; v < X.Length; v ++)
        {
			var i = v / size;
			var j = v - i * size;
			h[i, j] = X[v].y;
        }

		if (Input.GetKeyDown ("r")) 
		{
			//TODO: Add random water.
			var r = Random.Range(0, 1.0f);
			var i = Random.Range(0, size);
			var j = Random.Range(0, size);
			h[i, j] += r;
			// compute surround sum and deduce h with r / sum
			int sum = 0;
			int i_b = Mathf.Max(0, i - 1);
			int i_e = Mathf.Min(i + 1, size - 1);
			int j_b = Mathf.Max(0, j - 1);
			int j_e = Mathf.Min(j + 1, size - 1);
			for(int i_s = i_b; i_s <= i_e; i_s ++)
            {
				for(int j_s = j_b; j_s <= j_e; j_s ++)
                {
					if (i_s != i && j_s != j)
						sum ++;
                }
            }
			var dr = r / sum;
			for(int i_s = i_b; i_s <= i_e; i_s ++)
            {
				for(int j_s = j_b; j_s <= j_e; j_s ++)
                {
					if (i_s != i && j_s != j)
						h[i_s, j_s] -= dr;
                }
            }
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

        //TODO: Store h back into X.y and recalculate normal.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                var v = i * size + j;
                X[v].y = h[i, j];
            }
        }
        mesh.vertices = X;
		mesh.RecalculateNormals();
	}
}
