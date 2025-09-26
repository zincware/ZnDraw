import { createBrowserRouter, RouterProvider } from 'react-router-dom';
import RoomPage from './pages/roomPage';
import MainPage from './pages/landingPage';


const LandingPage = () => {
  return (
    <div className="app-layout">
      <header>
        <h1>Welcome to Collaborative Viewer</h1>
        <p>Please select or create a room to join.</p>
        <a href="/room/testroom">Join Test Room</a>
      </header>
    </div>
  );
};

const router = createBrowserRouter([
  {
    path: '/',
    element: <LandingPage />, // A landing page with instructions or a room list
  },
  {
    path: '/room/:roomId', // The main viewer route
    element: <RoomPage />,
  },
  {
    path: '/main',
    element: <MainPage />,
  },
]);

export function App() {
  return <RouterProvider router={router} />;
}