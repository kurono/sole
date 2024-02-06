// ES6

/* By Tsivilskiy Ilya, 10.08.2018 */

/* Class definitions */

/**
 * The Matrix
 */
class Matrix extends Array {
    constructor() {
        super();
    }

    get nRows() {
        return this.length;
    }

    get nCols() {
        let result = 0;
        if (this.length !== 0) {
            result = this[0].length;
        }
        return result;
    }

    get size() {
        return {"nRows": this.nRows, "nCols": this.nCols};
    }

    static multiply(a, b) {
        let ri, ci, i;

        let anr = a.nRows, anc = a.nCols,
            bnr = b.nRows, bnc = b.nCols,
            m = Matrix.zeros(bnr, bnc);
        for (ri = 0; ri <  anr; ri++) {
            for (ci = 0; ci < bnc; ci++) {
                for (i = 0; i < anc; i++) {
                    m[ri][ci] = m[ri][ci] + a[ri][i] * b[i][ci];
                }
            }
        }

        return m;
    }

    static hasTheSameDimensions(a, b) {
        return (JSON.stringify(a.size) === JSON.stringify(b.size));
    }

    static scaleBy(a, n) {
        let ri, ci;
        let nr = a.nRows, nc = a.nCols,
            m = Matrix.zeros(nr, nc);
        for (ri = 0; ri <  nr; ri++) {
            for (ci = 0; ci < nc; ci++) {
                m[ri][ci] = a[ri][ci] * n;
            }
        }
        return m;
    }

    static shiftBy(a, n) {
        let ri, ci;
        let nr = a.nRows, nc = a.nCols,
            m = Matrix.zeros(nr, nc);
        for (ri = 0; ri <  nr; ri++) {
            for (ci = 0; ci < nc; ci++) {
                m[ri][ci] = a[ri][ci] + n;
            }
        }
        return m;
    }

    static add(a, b) {
        if (!Matrix.hasTheSameDimensions(a, b)) throw "Matrix dimensions must agree!";

        let ri, ci;
        let nr = a.nRows, nc = a.nCols,
            m = Matrix.zeros(nr, nc);
        for (ri = 0; ri <  nr; ri++) {
            for (ci = 0; ci < nc; ci++) {
                m[ri][ci] = a[ri][ci] + b[ri][ci];
            }
        }
        return m;
    }

    static subtract(a, b) {
        if (!Matrix.hasTheSameDimensions(a, b)) throw "Matrix dimensions must agree!";

        let negb = Matrix.scaleBy(Matrix.clone(b), -1);
        return Matrix.add(a, negb);
    }

    static transpose(a) {
        let m = Matrix.zeros(a.nCols, a.nRows);
        let ri, ci;
        let nr = m.nRows, nc = m.nCols;
        for (ri = 0; ri < nr; ri++) {
            for (ci = 0; ci < nc; ci++) {
                m[ri][ci] = a[ci][ri];
            }
        }
        return m;
    }

    static zeros(n, m) {
        let a = new Matrix();
        let ri, ci;
        let row;
        for (ri = 0; ri < n; ri++) {
            row = [];
            for (ci = 0; ci < m; ci++) {
                row.push(0);
            }
            a.push(row);
        }
        return a;
    }

    static random(n, m) {
        let a = new Matrix();
        let ri, ci;
        let row;
        for (ri = 0; ri < n; ri++) {
            row = [];
            for (ci = 0; ci < m; ci++) {
                row.push(Math.random());
            }
            a.push(row);
        }
        return a;
    }

    static clone(a) {
        let m = new Matrix();
        let ri, ci;
        let nr = a.nRows, nc = a.nCols;
        let row;
        for (ri = 0; ri < nr; ri++) {
            row = [];
            for (ci = 0; ci < nc; ci++) {
                row.push(a[ri][ci]);
            }
            m.push(row);
        }
        return m;
    }
}

/**
 * System of linear equations
  */
class SOLE {
    static get GAUSS_ELIMINATION() {
        return "gauss";
    }

    static get JACOBI_ITERATION() {
        return "jacobi";
    }

    static get GAUSS_SEIDEL_RELAXATION() {
        return "gauss-seidel";
    }

    constructor(A, B) {
        this.A = A;
        this.B = B;
        this.X = Matrix.zeros(A.nCols, 1);
        this.nIter = 0;
    }

    /**
     *
     * @returns {A * X - B}
     */
    get checkSolutionAsVector() {
        let diff = Matrix.zeros(this.X.nRows, 1);
        let AX = Matrix.multiply(this.A, this.X);
        diff = Matrix.subtract(AX, this.B);

        return diff;
    }

    /**
     *
     * @returns {sum of 'A * X - B'}
     */
    get checkSolutionAsScalar() {
        let diff = this.checkSolutionAsVector;
        let mag = 0;
        let ri;
        let nr = diff.nRows;
        for (ri=0; ri < nr; ri++) {
            mag = mag + Math.pow(diff[ri][0], 2);
        }
        mag = Math.sqrt(mag);

        return mag;
    }

    solve(solverType) {
        switch (solverType){
            case SOLE.GAUSS_ELIMINATION:
                this.solveByGaussElimination();
                break;
            case SOLE.JACOBI_ITERATION:
                this.solveByJacobiIteration();
                break;
            case SOLE.GAUSS_SEIDEL_RELAXATION:
                this.solveByGaussSeidelRelaxation();
                break;
            default :
                break;
        }
    }

    /*
     * Direct solution by Gauss elimination
     */	
	solveByGaussElimination() {
        let n = this.X.length;

        let a = Matrix.clone(this.A);
        let b = Matrix.clone(this.B);
        let pr, max, i, j, temp, t, alpha, sum;
		let EPSILON = 1e-05;
		
		for (pr = 0; pr < n; pr++) {
            // find pivot row and swap
            max = pr;
            for (i = pr + 1; i < n; i++) {
                if (Math.abs(a[i][pr]) > Math.abs(a[max][pr])) {
                    max = i;
                }
            }
            temp = a[pr]; a[pr] = a[max]; a[max] = temp;
            t = b[pr]; b[pr] = b[max]; b[max] = t;

            // singular or nearly singular
            if (Math.abs(a[pr][pr]) <= EPSILON) {
                throw "Matrix is singular or nearly singular";
            }

            // pivot within a and b
            for (i = pr + 1; i < n; i++) {
                alpha = a[i][pr] / a[pr][pr];
                b[i] -= alpha * b[pr];
                for (j = pr; j < n; j++) {
                    a[i][j] -= alpha * a[pr][j];
                }
            }
        }

        // back substitution
        let x = Matrix.clone(this.X);
		x = Matrix.scaleBy(x, 0.0);
        for (i = n - 1; i >= 0; i--) {
            sum = 0.0;
            for (j = i + 1; j < n; j++) {
                sum += a[i][j] * x[j][0];
            }
            x[i][0] = (b[i] - sum) / a[i][i];
        }
		
		this.X = x;
		
		console.log(this.B.size, this.X.size);
    }

    solveByJacobiIteration() {
        //throw "The method is not implemented yet!";

        let a = Matrix.clone(this.A);
        let b = Matrix.clone(this.B);

        let prev_x = Matrix.clone(this.X);
        let current_x = Matrix.clone(this.X);

        let ri, ci;
        let nr = a.nRows, nc = a.nCols;

        let x;

        let i, n = this.nIter;

        for (i=0; i < n; i++) {

            for (ri=0; ri<nr; ri++) {
                x = b[ri][0];

                for (ci=0; ci<nc; ci++) {
                    if (ci !== ri) {
                        x = x - a[ri][ci]*prev_x[ci][0];
                    }
                }
                x = x / a[ri][ri];

                current_x[ri][0] = x;
            }

            prev_x = Matrix.clone(current_x);
        }

        this.X = Matrix.clone(current_x);
    }
}

/* Main program starts here */

const N_ITER = 10;
const INITIAL_GUESS = 0.001;

let sys;

function log(...args) {
    console.log(...args);
}

/**
 * Convert HTML table to A and B matrices
 * @param tableId
 * @returns {{A, B}}
 */
function tableToMatrices(tableId) {
    let table = document.getElementById(tableId);

    let A = new Matrix();
    let B = new Matrix();

    let row = [];
    let cc;
    let ri, ci;
    let nr = table.rows.length;
    let nc;

    for (ri = 0; ri < nr; ri++) {
        nc = table.rows[ri].cells.length;
        for (ci = 0; ci < nc; ci++) {
            cc = parseFloat(table.rows[ri].cells[ci].innerHTML);
            if (ci < nc-1) {
                row.push(cc);
            } else {
                A.push(row);
                B.push([cc]);
                row = [];
            }
        }
    }

    return {A, B};
}

function tableToInitialGuess(tableId) {
    let table = document.getElementById(tableId);

    let v = new Matrix();

    let row = [];
    let cc;
    let ri, ci;
    let nr = table.rows.length;
    let nc;

    for (ri = 0; ri < nr; ri++) {
        nc = table.rows[ri].cells.length;
        row = [];
        for (ci = 0; ci < nc; ci++) {
            cc = parseFloat(table.rows[ri].cells[ci].innerHTML);
            row.push(cc);
        }
        v.push(row);
    }

    return v;
}

function solutionToTable(matrix, tableId) {
    let table = document.getElementById(tableId);
    let ri, ci;
    let nr = table.rows.length;
    let nc;
    for (ri = 0; ri < nr; ri++) {
        nc = table.rows[ri].cells.length;
        for (ci = 0; ci < nc; ci++) {
            table.rows[ri].cells[ci].innerHTML = matrix[ri][ci].toString();
        }
    }
}

function calculate(solverType) {
    let mxs = tableToMatrices("AnB");
    console.log("A=", mxs.A, "B=", mxs.B);

    /*let randA = Matrix.random(3, 3);
    let randB = Matrix.random(3, 1);
    sys = new SOLE(randA, randB);*/

    sys = new SOLE(mxs.A, mxs.B);

    switch (solverType) {
        case SOLE.GAUSS_ELIMINATION:
            sys.nIter = 0;
            break;
        case SOLE.JACOBI_ITERATION:
            sys.nIter = N_ITER;
            sys.X = tableToInitialGuess("X");
            break;
    }

    log("initial X=", sys.X);

    sys.solve(solverType);

    log("solved X=", sys.X);

    solutionToTable(sys.X, "X");

    log("root check  = ", sys.checkSolutionAsVector);

    document.getElementById("diff").innerHTML = "|| AX - B || = " + sys.checkSolutionAsScalar.toString();
}

function onButtonClick(solverType) {
    calculate(solverType);
}

function onResetButtonClick() {
    if (sys === undefined) return;

    let n = sys.X.nRows;
    sys.X = Matrix.shiftBy(Matrix.zeros(n, 1), INITIAL_GUESS);

    solutionToTable(sys.X, "X");

    document.getElementById("diff").innerHTML = "|| AX - B || = " + sys.checkSolutionAsScalar.toString();
}