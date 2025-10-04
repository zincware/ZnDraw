import { useFigureList, useDeleteFigure } from '../hooks/useFigures';

export default function FigureList() {
  const { data, isLoading, error } = useFigureList();
  const deleteFigureMutation = useDeleteFigure();

  if (isLoading) return <div>Loading figures...</div>;
  if (error) return <div>Error: {error.message}</div>;

  const handleDelete = (key: string) => {
    if (window.confirm(`Are you sure you want to delete figure "${key}"?`)) {
        deleteFigureMutation.mutate(key);
    }
  }

  return (
    <aside>
      <h2>Figures</h2>
      <ul>
        {data?.figures.map((key) => (
          <li
            key={key}
            style={{
              fontWeight: 'bold',
              cursor: 'pointer',
            }}
          >
           {/* <span onClick={() => setSelectedFigureKey(key)}>{key}</span> */}
           <button onClick={() => handleDelete(key)} style={{marginLeft: '10px'}}>X</button>
          </li>
        ))}
      </ul>
    </aside>
  );
}