#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <stmesh/stmesh.hpp>
#include <string>

int factorial(int input) noexcept {
  int result = 1;

  while (input > 0) {
    result *= input;
    --input;
  }

  return result;
}

void test_vtk() {
  const vtkNew<vtkNamedColors> colors;

  // Parse command line arguments.
  const std::string filename = "test.vtu";

  const vtkNew<vtkPoints> points;
  points->InsertNextPoint(0, 0, 0);
  points->InsertNextPoint(1, 0, 0);
  points->InsertNextPoint(1, 1, 0);
  points->InsertNextPoint(0, 1, 1);

  const vtkNew<vtkTetra> tetra;

  tetra->GetPointIds()->SetId(0, 0);
  tetra->GetPointIds()->SetId(1, 1);
  tetra->GetPointIds()->SetId(2, 2);
  tetra->GetPointIds()->SetId(3, 3);

  const vtkNew<vtkCellArray> cellArray;
  cellArray->InsertNextCell(tetra);

  const vtkNew<vtkUnstructuredGrid> unstructuredGrid;
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(VTK_TETRA, cellArray);

  // Write file.
  const vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();

  // Read and display file for verification that it was written correclty.
  const vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();

  const vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputConnection(reader->GetOutputPort());

  const vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

  const vtkNew<vtkRenderer> renderer;
  const vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("WriteVTU");

  const vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("cobalt_green").GetData());

  renderWindow->Render();
  renderWindowInteractor->Start();
}
