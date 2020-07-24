import core.stdc.stdio: printf;
import core.stdc.float_: FLT_EPSILON, DBL_EPSILON;
import std.math: abs;


alias double DEFAULT_SCALAR;

template Matrix(int M, int N, T = DEFAULT_SCALAR)
{
	alias T[N][M] mat;
	
	//U is the type of the optional inverse return value
	void rowReduce(U = typeof(null))(scope const mat val, ref mat oo, U* inv = null)
	{
		const bool Tall = M > N;
		const size_t PrelimRowMax = Tall ? N : M;
		const size_t ExcessRows = Tall ? M-N : 0;
		alias T[PrelimRowMax][PrelimRowMax] sqMat;
		
		static if (N <= 0 || M <= 0)
		{
			printf("Matrix must have non-zero size!\n");
			return;
		}
		
		immutable bool calcInverse = M == N && is(U == mat);
		
		if (calcInverse && inv is null)
		{
			printf("Inverse pointer must be non-null!\n");
			return;
		}
	
		oo = val;
		size_t columnOffset = 0;
		static if (calcInverse)
		{
			foreach (size_t rr, T[N] row; *inv)
			{
				foreach (size_t cc, T ent; row)
				{
					(*inv)[rr][cc] = rr==cc ? 1.0 : 0.0;
				}
			}
		}
		
		T[N] swapSpace;
		
		void swapRow(size_t a, size_t b)
		{
			void action(ref mat subj)
			{
				swapSpace = subj[a];
				subj[a] = subj[b];
				subj[b] = swapSpace;
			}
			
			action(oo);
			static if (calcInverse)
				action(*inv);
		}
		
		void inPlaceMult(size_t a, T scalar)
		{
			oo[a][] *= scalar;
			static if (calcInverse)
				(*inv)[a][] *= scalar;
		}
		
		void rowAddMult(size_t targ, size_t source, T scalar)
		{
			oo[targ][] += scalar*oo[source][];
			static if (calcInverse)
				(*inv)[targ][] += scalar*(*inv)[source][];
		}
		
		//List of leading entries
		size_t[2][PrelimRowMax] leadingEntries;
		size_t lEntryLength = 0;
		
		void pushLeadingEntry(size_t r, size_t c)
		{
			leadingEntries[lEntryLength] = [r,c];
			lEntryLength++;
		}
		
		immutable T EPSILON = is(T == float) ? FLT_EPSILON : DBL_EPSILON;
		T maxElement = oo[0][0];
		
		foreach (T[N] row; oo)
		{
			foreach (T num; row)
			{
				if (num > maxElement)
					maxElement = num;
			}

		}
		
		//Tests for approximate equality
		bool eql(T a, T b, T tolerance = EPSILON*maxElement)
		{
			return abs(a - b) < tolerance;
		}
		
		for (size_t row = 0; row < M && row+columnOffset < N; row++)
		{
			//'Row-subset' [row..M] is the set of remaining rows that don't include the previously processed rows
			size_t col = row+columnOffset;
			
			//Determine if the column (row+offset) in the row-subset contains a non zero entry
			size_t unitrow = row;
			while (unitrow < M && eql(oo[unitrow][col], 0.0))
			{
				oo[unitrow][col] = 0.0;
				unitrow++;
			}
				
			if (unitrow < M)
			{
			//If so,
				//move it to the top of this row-subset and unit-reduce its leading entry
				swapRow(row, unitrow);
				inPlaceMult(row, 1.0/oo[row][col]);
				oo[row][col] = 1.0;
				//Store this coordinate set in a list of leading units
				pushLeadingEntry(row,col);
				//Eliminate each of the other non-zero column entries in [1+row..M]
				for (size_t remEnt = row+1; remEnt < M; remEnt++)
				{
					//Do this by reducing the subset-top row to have a leading unit entry
					//(manually set leading entry and anything it equals in that row to 0 or 1 respectivly)
					//Use that to zero out the rest of the column, when the entries aren't zero already.
					T entVal = oo[remEnt][col];
					if (!eql(entVal, 0.0))
						rowAddMult(remEnt,row,-entVal);
					
					oo[remEnt][col] = 0.0;
				}
			}
			else
			{
			//If the remaining column is zero-only,
				//increase the columnOffset by one, decrease row (to keep it the same) and continue
				columnOffset++;
				row--;
			}
		}
		
		//Eliminate all the 'tops' of the columns with leading unit entries
		for (size_t ii = 0; ii<lEntryLength; ii++)
		{
			size_t r = leadingEntries[ii][0];
			size_t c = leadingEntries[ii][1];
			
			for (size_t higherRow = 0; higherRow < r; higherRow++)
			{
				T entVal = oo[higherRow][c];
				if (!eql(entVal, 0.0))
					rowAddMult(higherRow, r, -entVal);
				oo[higherRow][c] = 0.0;
			}
		}
		
	}
	
	void print(mat a)
	{
		printf("[");
		
		foreach (T[N] row; a)
		{
			printf("[");
			foreach (T val; row)
			{
				printf("%f, ", val);
			}
			printf("]\n");
		}
		
		printf("]\n");
	}
}

//void MatrixMult(int M, int N, int P, T)(ref T[N][M] a, ref T[P][N] b, ref T[P][M] c)
void MatrixMult(int M, int N, int P, T = DEFAULT_SCALAR)(ref Matrix!(M,N,T).mat a, ref Matrix!(N,P,T).mat b, ref Matrix!(M,P,T).mat c)
{
	for (size_t rr = 0; rr < M; rr++)
	{
		for (size_t cc = 0; cc < P; cc++)
		{
			T acclum = 0.0;
			for (size_t ii = 0; ii < P; ii++)
				acclum += a[rr][ii]*b[ii][cc];
			
			c[rr][cc] = acclum;
		}
	}
}

alias Matrix!(3,3) matrix;

extern(C):
void main()
{
	
	matrix.mat a = 0, b = 0, c = 0, d = 0;
	a[0] = [-1.0, 0.0, 0.0];
	a[1] = [0.0, 2.0, 1.0];
	a[2] = [0.0, 0.0, 0.125];
	matrix.rowReduce(a,b,&c);
	printf("Matrix A (original):\n");
	matrix.print(a);
	printf("Matrix B (reduced):\n");
	matrix.print(b);
	printf("Matrix C (inverted):\n");
	matrix.print(c);
	MatrixMult!(3,3,3)(a,c,d);
	printf("Matrix D (multiplied):\n");
	matrix.print(d);
}