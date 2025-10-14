import React, { Component, ReactNode } from "react";

interface Props {
  geometryKey: string;
  children: ReactNode;
}

interface State {
  hasError: boolean;
  error: Error | null;
}

/**
 * Error boundary that wraps individual geometry components.
 * Prevents one geometry from crashing the entire scene.
 *
 * Following SOLID principles - isolates geometry failures to prevent cascade.
 */
export class GeometryErrorBoundary extends Component<Props, State> {
  constructor(props: Props) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error: Error): State {
    return { hasError: true, error };
  }

  componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
    console.error(
      `Error in geometry "${this.props.geometryKey}":`,
      error,
      "\nStack:",
      error.stack,
      "\nComponent Stack:",
      errorInfo.componentStack
    );
  }

  render() {
    if (this.state.hasError) {
      // Return null to render nothing instead of crashing the app
      // The geometry will simply not appear in the scene
      console.warn(
        `Geometry "${this.props.geometryKey}" failed to render and was removed from the scene.`
      );
      return null;
    }

    return this.props.children;
  }
}
