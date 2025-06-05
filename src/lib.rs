use numpy::{PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::{exceptions::PyValueError, prelude::*, BoundObject};
use rayon::prelude::*;

pub mod gms_nav;
pub mod goes_nav;

use gms_nav::GMSNavigation;
use goes_nav::GOESNavigation;

fn get_lonlat_impl<'py, F>(
    py: Python<'py>,
    lin: PyReadonlyArray2<'_, f64>,
    ele: PyReadonlyArray2<'_, f64>,
    f: F,
) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>)>
where
    F: Fn(f64, f64) -> (f64, f64) + Sync,
{
    let dim = lin.dims();
    let xdim = dim[0];
    let ydim = dim[1];
    let lin_slice = lin.as_slice()?;
    let ele_slice = ele.as_slice()?;
    let lon_py = unsafe { PyArray2::<f64>::new(py, [xdim, ydim], false) };
    let lat_py = unsafe { PyArray2::<f64>::new(py, [xdim, ydim], false) };
    let lon_buf = unsafe { lon_py.as_slice_mut()? };
    let lat_buf = unsafe { lat_py.as_slice_mut()? };

    py.allow_threads(|| {
        lon_buf
            .par_iter_mut()
            .zip(lat_buf.par_iter_mut())
            .enumerate()
            .for_each(|(k, (lon_slot, lat_slot))| {
                let lin_v = lin_slice[k];
                let ele_v = ele_slice[k];
                let (lat_v, lon_v) = f(lin_v, ele_v);
                *lon_slot = lon_v;
                *lat_slot = lat_v;
            });
    });

    Ok((lon_py.into_bound(), lat_py.into_bound()))
}

#[pyclass(name = "GOESNav", module = "satnav")]
pub struct PyGOESNavigation {
    inner: GOESNavigation,
}

#[pymethods]
impl PyGOESNavigation {
    #[new]
    #[pyo3(signature = (navblock))]
    pub fn py_new(navblock: PyReadonlyArray1<i32>) -> PyResult<Self> {
        let data = navblock.to_vec().unwrap();
        let inner = GOESNavigation::new(&data);
        Ok(PyGOESNavigation { inner })
    }

    pub fn get_subpoint(&self) -> (f64, f64) {
        self.inner.get_subpoint()
    }

    pub fn get_lonlat<'py>(
        &self,
        py: Python<'py>,
        lin: PyReadonlyArray2<'_, f64>,
        ele: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>)> {
        get_lonlat_impl(py, lin, ele, |lin, ele| {
            self.inner.image_to_geographic(lin, ele)
        })
    }

    pub fn geographic_to_image(&self, lon: f64, lat: f64) -> PyResult<(f64, f64)> {
        let (lin, ele) = self.inner.geographic_to_image(lon, lat);
        Ok((lin, ele))
    }
}

#[pyclass(name = "GMSNav", module = "satnav")]
pub struct PyGMSNavigation {
    inner: GMSNavigation,
}

#[pymethods]
impl PyGMSNavigation {
    #[new]
    #[pyo3(signature = (
        atit,
        orbt1,
        rstep,
        rsamp,
        rfcl,
        rfcp,
        dsct,
        dspin,
        sens,
        vmis,
        elmis,
    ))]
    fn py_new(
        atit: PyReadonlyArray2<f64>,
        orbt1: PyReadonlyArray2<f64>,
        rstep: f64,
        rsamp: f64,
        rfcl: f64,
        rfcp: f64,
        dsct: f64,
        dspin: f64,
        sens: f64,
        vmis: PyReadonlyArray1<f64>,
        elmis: PyReadonlyArray2<f64>,
    ) -> PyResult<Self> {
        /* ---------------- shape checks ---------------- */
        if atit.shape() != [10, 33] {
            return Err(PyValueError::new_err("ATIT must have shape (10, 33)"));
        }
        if orbt1.shape() != [35, 18] {
            return Err(PyValueError::new_err("ORBT1 must have shape (35, 18)"));
        }
        if vmis.len() != 3 {
            return Err(PyValueError::new_err("VMIS must have length 3"));
        }
        if elmis.shape() != [3, 3] {
            return Err(PyValueError::new_err("ELMIS must have shape (3, 3)"));
        }

        /* -------- copy NumPy buffers into fixed-size Rust arrays -------- */
        let mut atit_rs = [[0.0_f64; 33]; 10];
        for (i, row) in atit.as_array().outer_iter().enumerate() {
            atit_rs[i].copy_from_slice(row.as_slice().unwrap());
        }

        let mut orbt_rs = [[0.0_f64; 18]; 35];
        for (i, row) in orbt1.as_array().outer_iter().enumerate() {
            orbt_rs[i].copy_from_slice(row.as_slice().unwrap());
        }

        let vmis_rs = [
            vmis.get([0]).unwrap().clone(),
            vmis.get([1]).unwrap().clone(),
            vmis.get([2]).unwrap().clone(),
        ];

        let mut elmis_rs = [[0.0_f64; 3]; 3];
        for (i, row) in elmis.as_array().outer_iter().enumerate() {
            elmis_rs[i].copy_from_slice(row.as_slice().unwrap());
        }

        /* -------------- create inner Rust navigator -------------------- */
        let inner = GMSNavigation::new(
            atit_rs, orbt_rs, rstep, rsamp, rfcl, rfcp, dsct, dspin, sens, vmis_rs, elmis_rs,
        );
        Ok(PyGMSNavigation { inner })
    }

    pub fn get_lonlat<'py>(
        &self,
        py: Python<'py>,
        lin: PyReadonlyArray2<'_, f64>,
        ele: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>)> {
        get_lonlat_impl(py, lin, ele, |lin, ele| {
            self.inner.image_to_geographic(lin, ele)
        })
    }

    pub fn geographic_to_image(&self, lon: f64, lat: f64) -> PyResult<(f64, f64)> {
        let (lin, ele) = self.inner.geographic_to_image(lon, lat);
        Ok((lin, ele))
    }
}

/// Module definition ----------------------------------------------------------
#[pymodule]
fn satnav(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGOESNavigation>()?;
    m.add_class::<PyGMSNavigation>()?;
    Ok(())
}
